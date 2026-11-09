"""Generation-bound inputs for reproducible loading and resumable caches."""

import base64
import hashlib
import json
import os
import re
from collections.abc import Mapping
from functools import wraps
from http import HTTPStatus
from inspect import isgeneratorfunction
from urllib.parse import urlsplit

CRC32C_BYTES = 4


def normalize_identity(callset_path, identity):
    if isinstance(identity, Mapping) and not identity:
        return {}
    if not isinstance(identity, Mapping) or set(identity) != {'vcf', 'index'}:
        msg = 'Expected VCF and index identities'
        raise ValueError(msg)
    if not callset_path.endswith(('.vcf', '.vcf.gz', '.vcf.bgz')):
        msg = 'Bound callset must be a VCF'
        raise ValueError(msg)
    return {
        role: (
            None
            if role == 'index' and identity[role] is None
            else _normalize_object(identity[role], expected_path)
        )
        for role, expected_path in [
            ('vcf', callset_path),
            ('index', callset_path + '.tbi'),
        ]
    }


def _normalize_object(value, expected_path):
    required = {'path', 'generation', 'size', 'crc32c'}
    if (
        not isinstance(value, Mapping)
        or not required <= set(value)
        or set(value) - required - {'sha256'}
    ):
        msg = 'Malformed object identity'
        raise ValueError(msg)
    value = dict(value)
    if not isinstance(value['path'], str):
        msg = 'Invalid object path'
        raise ValueError(msg)  # noqa: TRY004 - Pydantic validation error
    uri = urlsplit(value['path'])
    if (
        value['path'] != expected_path
        or uri.scheme != 'gs'
        or not uri.netloc
        or not uri.path.strip('/')
        or uri.query
        or uri.fragment
        or '*' in value['path']
    ):
        msg = 'Identity must name the exact GCS VCF/index pair'
        raise ValueError(msg)
    for field in ('generation', 'size'):
        if not isinstance(value[field], str) or not re.fullmatch(
            r'[1-9][0-9]*',
            value[field],
        ):
            msg = 'Invalid generation or size'
            raise ValueError(msg)
    try:
        crc = base64.b64decode(value['crc32c'], validate=True)
    except (ValueError, TypeError) as exc:
        msg = 'Invalid CRC32C'
        raise ValueError(msg) from exc
    if len(crc) != CRC32C_BYTES or base64.b64encode(crc).decode() != value['crc32c']:
        msg = 'Invalid CRC32C'
        raise ValueError(msg)
    if 'sha256' in value and (
        not isinstance(value['sha256'], str)
        or not re.fullmatch(r'[0-9a-f]{64}', value['sha256'])
    ):
        msg = 'Invalid SHA256'
        raise ValueError(msg)
    return value


def verify_callset_identity(callset_path, identity, client=None):
    identity = normalize_identity(callset_path, identity)
    if not identity:
        reject_unbound_managed_alias(callset_path, client)
        return
    if client is None:
        from googleapiclient.discovery import build

        client = build('storage', 'v1', cache_discovery=False)
    for expected in identity.values():
        if expected is None:
            continue
        uri = urlsplit(expected['path'])
        # Inspect the live generation, not historical bytes. Existing pipeline
        # google-api-python-client uses application default credentials.
        actual = (
            client.objects()
            .get(
                bucket=uri.netloc,
                object=uri.path.lstrip('/'),
                ifGenerationMatch=expected['generation'],
            )
            .execute()
        )
        if any(
            str(actual.get(key)) != expected[key]
            for key in ('generation', 'size', 'crc32c')
        ):
            msg = 'Callset object changed after validation'
            raise ValueError(msg)


def identity_cache_key(callset_path, identity):
    normalized = normalize_identity(callset_path, identity)
    return hashlib.sha256(
        json.dumps(normalized, sort_keys=True, separators=(',', ':')).encode(),
    ).hexdigest()


def guard_callset(task_class):
    """Luigi's inherits copies parameters only, so wrap each loading task explicitly."""

    def verify(task):
        verify_callset_identity(task.callset_path, task.callset_identity)
        verify_run_binding(task)

    complete = task_class.complete

    @wraps(complete)
    def guarded_complete(self, *args, **kwargs):
        verify(self)
        result = complete(self, *args, **kwargs)
        verify(self)
        return result

    task_class.complete = guarded_complete

    run = task_class.run
    if isgeneratorfunction(run):

        @wraps(run)
        def guarded_run(self, *args, **kwargs):
            verify(self)
            generator = run(self, *args, **kwargs)
            sent = None
            while True:
                verify(self)
                try:
                    dependency = generator.send(sent)
                except StopIteration:
                    break
                sent = yield dependency
                verify(self)
            verify(self)
    else:

        @wraps(run)
        def guarded_run(self, *args, **kwargs):
            verify(self)
            result = run(self, *args, **kwargs)
            # BaseWriteTask.run returns only after Hail's write has materialized.
            verify(self)
            return result

    task_class.run = guarded_run
    return task_class


def _status(exc):
    return getattr(getattr(exc, 'resp', None), 'status', None)


def run_binding(task):
    def scalar(value):
        return getattr(value, 'value', value)

    return {
        'schema_version': 1,
        'callset_identity': normalize_identity(
            task.callset_path,
            task.callset_identity,
        ),
        'callset_path': task.callset_path,
        'reference_genome': scalar(task.reference_genome),
        'dataset_type': scalar(task.dataset_type),
        'sample_type': scalar(task.sample_type),
        'project_guids': sorted(task.project_guids),
        'validations_to_skip': sorted(task.validations_to_skip),
        'skip_check_sex_and_relatedness': task.skip_check_sex_and_relatedness,
        'skip_expect_tdr_metrics': task.skip_expect_tdr_metrics,
    }


def ensure_run_binding(run_path, binding, client=None):
    """Create once; never attach new provenance to existing unbound outputs."""
    payload = json.dumps(binding, sort_keys=True, separators=(',', ':')).encode()
    digest = hashlib.sha256(payload).hexdigest()
    manifest = run_path.rstrip('/') + '/_CALLSET_BINDING.json'
    mismatch = 'Run ID is already bound to different inputs or parameters'
    unbound = 'Existing unbound run outputs cannot be resumed with a new binding'
    if not manifest.startswith('gs://'):
        return _ensure_local_run_binding(manifest, payload, mismatch, unbound)
    return _ensure_gcs_run_binding(manifest, payload, digest, mismatch, unbound, client)


def _ensure_local_run_binding(manifest, payload, mismatch, unbound):
    import tempfile
    from pathlib import Path

    path = Path(manifest)
    if path.exists():
        if path.read_bytes() != payload:
            raise ValueError(mismatch)
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    if any(not p.name.startswith('.binding-') for p in path.parent.iterdir()):
        if path.exists() and path.read_bytes() == payload:
            return
        raise ValueError(unbound)
    # link() publishes a fully written file atomically without replacement.
    with tempfile.NamedTemporaryFile(dir=path.parent, prefix='.binding-') as temp:
        temp.write(payload)
        temp.flush()
        os.fsync(temp.fileno())
        try:
            os.link(temp.name, path)
        except FileExistsError:
            if path.read_bytes() != payload:
                raise ValueError(mismatch) from None
    return


def _matches_run_binding(objects, bucket, name, digest, mismatch):
    try:
        existing = objects.get(bucket=bucket, object=name).execute()
    except Exception as exc:
        if _status(exc) == HTTPStatus.NOT_FOUND:
            return False
        raise
    if existing.get('metadata', {}).get('binding_sha256') != digest:
        raise ValueError(mismatch)
    return True


def _ensure_gcs_run_binding(manifest, payload, digest, mismatch, unbound, client):
    if client is None:
        from googleapiclient.discovery import build

        client = build('storage', 'v1', cache_discovery=False)
    uri = urlsplit(manifest)
    bucket, name = uri.netloc, uri.path.lstrip('/')
    objects = client.objects()

    if _matches_run_binding(objects, bucket, name, digest, mismatch):
        return
    existing = objects.list(
        bucket=bucket,
        prefix=name.rsplit('/', 1)[0] + '/',
        maxResults=1,
    ).execute()
    if existing.get('items'):
        if _matches_run_binding(objects, bucket, name, digest, mismatch):
            return
        raise ValueError(unbound)
    from googleapiclient.http import MediaInMemoryUpload

    try:
        objects.insert(
            bucket=bucket,
            name=name,
            ifGenerationMatch='0',
            body={'metadata': {'binding_sha256': digest}},
            media_body=MediaInMemoryUpload(payload, mimetype='application/json'),
        ).execute()
    except Exception:
        # A timeout may have committed the write. Reconcile before any retry.
        if not _matches_run_binding(objects, bucket, name, digest, mismatch):
            raise


def verify_run_binding(task):
    if not task.callset_identity or not hasattr(task, 'run_id'):
        return
    from loading_pipeline.lib.paths import runs_path

    ensure_run_binding(
        os.path.join(runs_path(task.reference_genome, task.dataset_type), task.run_id),
        run_binding(task),
    )


def reject_unbound_managed_alias(callset_path, client=None):
    if os.environ.get('CALLSET_BINDING_ENABLED', '').lower() not in {'1', 'true'}:
        return
    if not callset_path.startswith('gs://') or not callset_path.endswith(
        ('.vcf', '.vcf.gz', '.vcf.bgz'),
    ):
        return
    if '*' in callset_path:
        import hailtop.fs as hfs

        paths = [entry.path for entry in hfs.ls(callset_path)]
    else:
        paths = [callset_path]
    if client is None:
        from googleapiclient.discovery import build

        client = build('storage', 'v1', cache_discovery=False)
    for path in paths:
        uri = urlsplit(path)
        actual = (
            client.objects()
            .get(bucket=uri.netloc, object=uri.path.lstrip('/'))
            .execute()
        )
        if 'seqr_binding' in (actual.get('metadata') or {}):
            msg = 'Managed callset aliases require a bound input identity'
            raise ValueError(msg)
