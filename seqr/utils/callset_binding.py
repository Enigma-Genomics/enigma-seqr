"""Bind validation to a concrete GCS object version, never to a mutable alias.

Opt-in for installations whose loading directory contains managed aliases. This
module snapshots legacy root objects with generation-pinned, create-only writes.
Managed aliases reuse their unique prepared objects. The unchanged loader receives
only that unique URI; it does not enforce the signed identity itself.
"""
import json
import fnmatch
import re
import uuid
from urllib.parse import urlsplit

from django.conf import settings
from django.core import signing
from google.cloud import storage

SALT = 'seqr.callset-validation.v1'
ERROR = 'Callset identity is missing, invalid, or changed. Re-stage this file if needed, then return to the VCF step and validate again.'


def _fail():
    raise ValueError(ERROR)


def _path(path):
    if not isinstance(path, str):
        _fail()
    parsed = urlsplit(path)
    if (parsed.scheme != 'gs' or not parsed.netloc or parsed.query or parsed.fragment
            or any(p in ('', '.', '..') for p in parsed.path.lstrip('/').split('/'))
            or re.search(r'[^A-Za-z0-9_./:-]', path)):
        _fail()
    return parsed.netloc, parsed.path.lstrip('/')


def enabled_for(path, client=None):
    if not getattr(settings, 'CALLSET_BINDING_ENABLED', False):
        return False
    root = settings.LOADING_DATASETS_DIR.rstrip('/') + '/'
    if not path.startswith(root):
        return False
    # Unmanaged sharded callsets and MatrixTables retain their existing flow.
    if '*' in path:
        # A broad legacy shard glob must not accidentally load a managed alias.
        relative_pattern = path[len(root):]
        if '/' not in relative_pattern or '*' in relative_pattern.split('/', 1)[0]:
            raise ValueError('Root VCF wildcards may include future mutable aliases. Select an individual file or a nested shard directory.')
        bucket, pattern = _path(path.replace('*', 'WILDCARD'))
        pattern = pattern.replace('WILDCARD', '*')
        for blob in (client or storage.Client()).list_blobs(bucket, prefix=pattern.split('*', 1)[0]):
            if not fnmatch.fnmatchcase(blob.name, pattern):
                continue
            if 'seqr_binding' in (blob.metadata or {}):
                _fail()
        return False
    if path.endswith('.mt'):
        return False
    _path(path)
    return path.endswith(('.vcf', '.vcf.gz', '.vcf.bgz'))


def _blob(client, path):
    bucket, name = _path(path)
    blob = client.bucket(bucket).get_blob(name)
    if blob is None:
        _fail()
    return blob


def _identity(path, blob):
    return {'path': path, 'generation': str(blob.generation), 'crc32c': blob.crc32c, 'size': str(blob.size)}


def verify_identity(identity, client=None):
    client = client or storage.Client()
    if not isinstance(identity, dict) or set(identity) != {'vcf', 'index'}:
        _fail()
    vcf = identity['vcf']
    index = identity['index']
    if (not isinstance(vcf, dict) or not isinstance(vcf.get('path'), str)
            or (index is not None and (not isinstance(index, dict) or index.get('path') != vcf['path'] + '.tbi'))):
        _fail()
    for expected in (vcf, index) if index is not None else (vcf,):
        if (not isinstance(expected.get('size'), str) or not re.fullmatch(r'[1-9][0-9]*', expected['size'])
                or not isinstance(expected.get('crc32c'), str) or not re.fullmatch(r'[A-Za-z0-9+/]{6}==', expected['crc32c'])):
            _fail()
        if not isinstance(expected.get('generation'), str) or not re.fullmatch(r'[1-9][0-9]*', expected['generation']):
            _fail()
        actual = _identity(expected.get('path'), _blob(client, expected.get('path')))
        if any(expected.get(k) != v for k, v in actual.items()):
            _fail()
    return identity


def _snapshot_legacy_root(identity, client):
    """Copy a selected version once; reconcile lost replies read-only, never overwrite."""
    root = settings.LOADING_DATASETS_DIR.rstrip('/') + '/'
    prefix = getattr(settings, 'CALLSET_SNAPSHOT_PREFIX', '.seqr-snapshots').strip('/')
    suffix = next(ext for ext in ('.vcf.gz', '.vcf.bgz', '.vcf') if identity['vcf']['path'].endswith(ext))
    destination = root + prefix + '/' + str(uuid.uuid4()) + '/input' + suffix
    _path(destination)
    snapshot = {'vcf': None, 'index': None}
    for kind in ('vcf', 'index'):
        original = identity[kind]
        if original is None:
            continue
        target_path = destination + ('.tbi' if kind == 'index' else '')
        bucket_name, source_name = _path(original['path'])
        bucket = client.bucket(bucket_name)
        source = bucket.blob(source_name, generation=int(original['generation']))
        _, target_name = _path(target_path)
        target = bucket.blob(target_name)
        marker = json.dumps(original, sort_keys=True, separators=(',', ':'))
        target.metadata = {'seqr_snapshot_source': marker}
        token = None
        try:
            while True:
                token, _, _ = target.rewrite(source, token=token, if_generation_match=0,
                    if_source_generation_match=int(original['generation']), retry=None)
                if token is None:
                    break
        except Exception:
            # A response can be lost after the write committed. Do not repeat an
            # uncertain operation; a matching completed unique object is enough.
            existing = bucket.get_blob(target_name)
            if existing is None:
                raise ValueError('Could not finish the callset snapshot. Return to the VCF step and validate again.') from None
            target = existing
        # Confirm actual persisted metadata, even after a successful rewrite.
        target = bucket.get_blob(target_name)
        if (target is None or (target.metadata or {}).get('seqr_snapshot_source') != marker
                or str(target.size) != original['size'] or target.crc32c != original['crc32c']):
            _fail()
        snapshot[kind] = _identity(target_path, target)
    verify_identity(snapshot, client)
    return snapshot


def resolve_callset(path, client=None):
    """Resolve an atomic alias binding, or pin a legacy nested unique input."""
    if not enabled_for(path, client):
        return path, None
    client = client or storage.Client()
    root = settings.LOADING_DATASETS_DIR.rstrip('/') + '/'
    blob = _blob(client, path)
    metadata = blob.metadata or {}
    raw_binding = metadata.get('seqr_binding')
    relative = path[len(root):]
    if 'seqr_binding' in metadata:
        try:
            binding = json.loads(raw_binding)
        except (TypeError, ValueError):
            _fail()
        if not isinstance(binding, dict) or binding.get('schema_version') != 1:
            _fail()
        accession = binding.get('accession_id')
        if not isinstance(accession, str) or not re.fullmatch(r'[A-Za-z0-9_-]+', accession):
            _fail()
        if relative != f'{accession}.gatk.vcf.gz':
            _fail()
        identity = {k: binding.get(k) for k in ('vcf', 'index')}
        vcf = identity['vcf']
        if not isinstance(vcf, dict):
            _fail()
        prefix = getattr(settings, 'CALLSET_PREPARED_PREFIX', 'lossless-v1').strip('/')
        expected_prefix = root + prefix + '/' + accession + '/'
        resolved = vcf.get('path')
        if not isinstance(resolved, str) or not resolved.startswith(expected_prefix):
            _fail()
        tail = resolved[len(expected_prefix):]
        if not re.fullmatch(r'[1-9][0-9]*-[1-9][0-9]*/[0-9a-fA-F-]{36}\.gatk\.vcf\.gz', tail):
            _fail()
        for expected in identity.values():
            if (not isinstance(expected, dict) or not isinstance(expected.get('sha256'), str)
                    or not re.fullmatch(r'[a-f0-9]{64}', expected['sha256'])):
                _fail()
        # The alias bytes and metadata must identify the same verified payload.
        if str(blob.size) != vcf.get('size') or blob.crc32c != vcf.get('crc32c'):
            _fail()
    else:
        # Legacy single files, including older root filenames, still receive
        # native generation binding. Present invalid alias metadata never falls back.
        resolved = path
        bucket, name = _path(path + '.tbi')
        index_blob = client.bucket(bucket).get_blob(name)
        identity = {'vcf': _identity(path, blob), 'index': _identity(path + '.tbi', index_blob) if index_blob else None}
    verify_identity(identity, client)
    if 'seqr_binding' not in metadata and '/' not in relative:
        identity = _snapshot_legacy_root(identity, client)
        resolved = identity['vcf']['path']
    return resolved, identity


def _context(request_json, user, path):
    return {
        'source': path, 'user': str(user.pk),
        'genomeVersion': request_json.get('genomeVersion'),
        'sampleType': request_json.get('sampleType'),
        'datasetType': request_json.get('datasetType', 'SNV_INDEL'),
    }


def sign_validation(request_json, user, path, identity, samples):
    return signing.dumps({
        **_context(request_json, user, path), 'identity': identity,
        'samples': sorted(samples),
    }, salt=SALT, compress=True)


def validated_callset(request_json, user, path, client=None):
    if not enabled_for(path, client):
        return path, None
    try:
        payload = signing.loads(request_json.get('callsetValidation', ''), salt=SALT,
                                max_age=getattr(settings, 'CALLSET_VALIDATION_MAX_AGE', 3600))
    except (signing.BadSignature, TypeError, ValueError):
        _fail()
    if not isinstance(payload, dict) or any(payload.get(k) != v for k, v in _context(request_json, user, path).items()):
        _fail()
    samples = request_json.get('vcfSamples')
    if not isinstance(samples, list) or not all(isinstance(sample, str) for sample in samples):
        _fail()
    if payload.get('samples') != sorted(samples):
        _fail()
    identity = verify_identity(payload.get('identity'), client)
    return identity['vcf']['path'], identity
