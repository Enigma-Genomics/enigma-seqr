# ruff: noqa: PLR2004
# Copyright 2026 Enigma Genomics
"""Offline intake of an attested joint trio VCF; never clinical preparation."""

import argparse
import gzip
import hashlib
import json
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path

import pysam


def require(ok, message):
    if not ok:
        raise ValueError(message)


def identifier(value):
    return isinstance(value, str) and re.fullmatch(
        r'[A-Za-z0-9][A-Za-z0-9_.-]{0,127}',
        value,
    )


def read_json(path):
    with Path(path).open('rb') as stream:
        data = stream.read(1_048_577)
    require(len(data) <= 1_048_576, 'JSON exceeds 1 MiB')
    return json.loads(data)


def contract(manifest, reference, family_id, joint_report_id):
    require(
        isinstance(manifest, dict) and manifest.get('schema_version') == 2,
        'Invalid intake schema',
    )
    require(
        manifest.get('kind') == 'attested_joint_trio_identity',
        'Not an attested joint trio',
    )
    for key, expected in [
        ('family_id', family_id),
        ('joint_report_id', joint_report_id),
    ]:
        require(
            identifier(expected) and manifest.get(key) == expected,
            f'{key} mismatch',
        )
    require(
        family_id != joint_report_id,
        'Family and report identities must be distinct',
    )
    require(
        set(manifest)
        == {
            'schema_version',
            'kind',
            'family_id',
            'joint_report_id',
            'genome_version',
            'sample_type',
            'dataset_type',
            'source_sha256',
            'header_sha256',
            'individuals',
            'joint_genotyping_attestation',
        },
        'Unexpected normalization manifest fields',
    )
    require(
        isinstance(manifest.get('header_sha256'), str)
        and re.fullmatch('[a-f0-9]{64}', manifest['header_sha256']),
        'Header hash required',
    )
    require(manifest.get('genome_version') == '38', 'Only GRCh38 supported')
    require(manifest.get('sample_type') in ('WES', 'WGS'), 'Invalid sample type')
    require(manifest.get('dataset_type') == 'SNV_INDEL', 'Invalid dataset type')
    require(
        re.fullmatch(r'[a-f0-9]{64}', str(manifest.get('source_sha256', ''))),
        'Source hash required',
    )
    provenance = manifest.get('joint_genotyping_attestation')
    require(
        isinstance(provenance, dict)
        and set(provenance) == {'confirmed', 'authority', 'reference', 'confirmed_by'},
        'Joint-genotyping attestation required',
    )
    require(provenance.get('confirmed') is True, 'Explicit confirmation required')
    require(
        provenance.get('authority') in ('user_confirmed', 'trusted_workflow'),
        'Invalid attestation authority',
    )
    for key in ('reference', 'confirmed_by'):
        require(
            isinstance(provenance.get(key), str)
            and 0 < len(provenance[key]) <= 500
            and provenance[key].strip() == provenance[key],
            f'Invalid provenance {key}',
        )
    require(
        isinstance(reference, dict) and reference.get('genome_version') == '38',
        'Reference build mismatch',
    )
    require(
        reference.get('verified') is True
        and isinstance(reference.get('verification_reference'), str)
        and bool(reference['verification_reference'].strip()),
        'Separately verified reference dictionary required',
    )
    contigs = reference.get('contigs')
    require(
        isinstance(contigs, dict) and 0 < len(contigs) <= 10000,
        'Invalid reference dictionary',
    )
    require(
        all(
            isinstance(k, str)
            # Broad GRCh38 includes HLA contigs such as HLA-A*01:01:01:01.
            and re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.:*-]{0,127}', k)
            and type(v) is int
            and 0 < v < 2**29
            for k, v in contigs.items()
        ),
        'Invalid reference contig',
    )
    require(
        contigs.get('1', contigs.get('chr1')) == 248956422,
        'GRCh38 chromosome 1 required',
    )
    return identity_mapping(manifest), contigs


def identity_mapping(manifest):
    people = manifest.get('individuals')
    require(
        isinstance(people, list)
        and len(people) == 3
        and all(isinstance(p, dict) for p in people),
        'Exactly three individuals required',
    )
    require(
        {p.get('role') for p in people} == {'proband', 'mother', 'father'},
        'Explicit roles required',
    )
    for field in ('sample_id', 'accession_id', 'source_sample_id'):
        values = [p.get(field) for p in people]
        require(
            all(identifier(v) for v in values) and len(set(values)) == 3,
            f'Invalid or duplicate {field}',
        )
    by_role = {p['role']: p for p in people}
    for person in people:
        require(
            set(person)
            == {
                'role',
                'sample_id',
                'accession_id',
                'source_sample_id',
                'mother_id',
                'father_id',
            },
            'Normalization accepts identity fields only; clinical metadata belongs in preparation review',
        )
        for parent in ('mother', 'father'):
            expected = (
                by_role[parent]['sample_id'] if person['role'] == 'proband' else None
            )
            require(
                parent + '_id' in person and person[parent + '_id'] == expected,
                'Pedigree identity mismatch',
            )
    return {p['source_sample_id']: p['sample_id'] for p in people}


def manifest_from_request(request, source_sha256):
    """Bridge an immutable identity-only queue to local intake, without clinical defaults.

    Caller downloads the exact pinned generation, hashes it, and separately supplies
    the operator-configured verified reference dictionary to normalize().
    """
    keys = {
        'schema_version',
        'kind',
        'intake_request_id',
        'source',
        'header_sha256',
        'joint_report_id',
        'family_id',
        'subject_key',
        'genome_version',
        'sample_type',
        'dataset_type',
        'individuals',
        'joint_genotyping_attestation',
    }
    require(
        isinstance(request, dict) and set(request) == keys,
        'Unexpected queue fields',
    )
    require(
        request['schema_version'] == 1 and request['kind'] == 'normalize_joint_trio',
        'Unsupported queue schema',
    )
    content = {
        key: value for key, value in request.items() if key != 'intake_request_id'
    }
    digest = hashlib.sha256(
        json.dumps(
            content,
            sort_keys=True,
            separators=(',', ':'),
            ensure_ascii=False,
        ).encode('utf-8'),
    ).hexdigest()
    require(
        request['intake_request_id'] == digest,
        'Intake request fingerprint mismatch',
    )
    require(identifier(request['joint_report_id']), 'Invalid joint report identifier')
    identity = hashlib.sha256(request['joint_report_id'].encode('utf-8')).hexdigest()
    require(
        request['family_id'] == 'JR_' + identity
        and request['subject_key'] == 'family:' + identity,
        'Derived family namespace mismatch',
    )
    source = request['source']
    require(
        isinstance(source, dict) and set(source) == {'bucket', 'path', 'generation'},
        'Invalid source pin',
    )
    for key, bound in (('bucket', 222), ('path', 1024), ('generation', 30)):
        require(
            isinstance(source[key], str)
            and 0 < len(source[key]) <= bound
            and not any(ord(c) < 32 for c in source[key]),
            'Invalid source pin',
        )
    require(
        re.fullmatch('[1-9][0-9]*', source['generation']),
        'Invalid pinned generation',
    )
    require(
        isinstance(source_sha256, str) and re.fullmatch('[a-f0-9]{64}', source_sha256),
        'Downloaded source hash required',
    )
    identity_mapping(request)
    return {
        'schema_version': 2,
        'kind': 'attested_joint_trio_identity',
        **{
            key: request[key]
            for key in (
                'family_id',
                'joint_report_id',
                'genome_version',
                'sample_type',
                'dataset_type',
                'header_sha256',
                'individuals',
                'joint_genotyping_attestation',
            )
        },
        'source_sha256': source_sha256,
    }


def normalize(  # noqa: C901, PLR0912, PLR0915, PLR0913 - explicit streaming checks
    source,
    output_dir,
    manifest,
    reference,
    family_id,
    joint_report_id,
    *,
    max_source_bytes=536_870_912,
    max_expanded_bytes=536_870_912,
    max_line_bytes=4_194_304,
    max_header_bytes=8_388_608,
    max_output_bytes=536_870_912,
    max_scratch_bytes=1_073_741_824,
    timeout=1800,
):
    """Write a new local directory; remove all newly created files on failure."""
    limits = (
        max_source_bytes,
        max_expanded_bytes,
        max_line_bytes,
        max_header_bytes,
        max_output_bytes,
        max_scratch_bytes,
        timeout,
    )
    require(
        all(type(v) is int and v > 0 for v in limits),
        'Limits must be positive integers',
    )
    mapping, contigs = contract(manifest, reference, family_id, joint_report_id)
    source, output_dir = Path(source), Path(output_dir)
    require(
        source.is_file() and source.stat().st_size <= max_source_bytes,
        'Source missing or too large',
    )
    started = time.monotonic()

    def check():
        require(time.monotonic() - started < timeout, 'Walltime limit exceeded')
        if output_dir.exists():
            sizes = [p.stat().st_size for p in output_dir.iterdir()]
            require(sum(sizes) <= max_scratch_bytes, 'Scratch limit exceeded')
            require(all(s <= max_output_bytes for s in sizes), 'Output limit exceeded')

    def digest(path):
        result = hashlib.sha256()
        total = 0
        with path.open('rb') as stream:
            for chunk in iter(lambda: stream.read(1_048_576), b''):
                check()
                total += len(chunk)
                require(
                    total <= (max_source_bytes if path == source else max_output_bytes),
                    'Hash input byte limit exceeded',
                )
                result.update(chunk)
        return result.hexdigest()

    require(digest(source) == manifest['source_sha256'], 'Source hash mismatch')
    output_dir.mkdir(parents=False, exist_ok=False)
    try:
        output = output_dir / 'callset.vcf.gz'
        with source.open('rb') as raw:
            compressed = raw.read(2) == b'\x1f\x8b'
        opener = gzip.open if compressed else open
        expanded = header_bytes = records = 0
        seen_header = False
        declared = {}
        reference_seen = False
        last = (-1, 0)
        row_digest = hashlib.sha256()
        header_digest = hashlib.sha256()
        source_samples = []
        ranks = {}
        with opener(source, 'rb') as src, pysam.BGZFile(str(output), 'w') as dst:
            while True:
                line = src.readline(max_line_bytes + 1)
                if not line:
                    break
                expanded += len(line)
                check()
                require(
                    len(line) <= max_line_bytes and expanded <= max_expanded_bytes,
                    'Line or expanded byte limit exceeded',
                )
                # Conservative uncompressed reservation bounds BGZF and index growth.
                require(
                    expanded * 4 + 1_048_576
                    <= min(max_output_bytes, max_scratch_bytes),
                    'Output/scratch reservation exceeded',
                )
                text = line.decode('utf-8').rstrip('\r\n')
                require(text and '\x00' not in text, 'Empty or NUL-containing line')
                if expanded == len(line):
                    require(
                        text in ('##fileformat=VCFv4.2', '##fileformat=VCFv4.3'),
                        'VCF 4.2/4.3 first line required',
                    )
                if text.startswith('#'):
                    require(not seen_header, 'Header after sample header')
                    header_bytes += len(line)
                    require(
                        line.endswith(b'\n'),
                        'Header lines require a final newline',
                    )
                    header_digest.update(line)
                    require(header_bytes <= max_header_bytes, 'Header limit exceeded')
                    if text.startswith('##reference='):
                        require(not reference_seen, 'Duplicate reference')
                        require(
                            text.split('=', 1)[1].rsplit('/', 1)[-1]
                            in (
                                'GRCh38',
                                'GRCh38.fa',
                                'GRCh38.fasta',
                                'hg38',
                                'hg38.fa',
                                'hg38.fasta',
                                'Homo_sapiens_assembly38.fasta',
                            ),
                            'Reference build mismatch',
                        )
                        reference_seen = True
                    if text.startswith('##contig='):
                        match = re.fullmatch(
                            r'##contig=<ID=([^,>]+),length=(\d+)(?:,[^>]*)?>',
                            text,
                        )
                        require(match is not None, 'Unsupported contig declaration')
                        name, length = match.groups()
                        require(
                            name not in declared and contigs.get(name) == int(length),
                            'Contig mismatch/duplicate',
                        )
                        declared[name] = int(length)
                    if text.startswith('#CHROM'):
                        require(
                            header_digest.hexdigest() == manifest['header_sha256'],
                            'Pinned header hash mismatch',
                        )
                        columns = text.split('\t')
                        require(
                            columns[:9]
                            == [
                                '#CHROM',
                                'POS',
                                'ID',
                                'REF',
                                'ALT',
                                'QUAL',
                                'FILTER',
                                'INFO',
                                'FORMAT',
                            ],
                            'Invalid columns',
                        )
                        source_samples = columns[9:]
                        require(
                            len(source_samples) == 3
                            and set(source_samples) == set(mapping),
                            'Source sample set mismatch',
                        )
                        require(
                            declared == contigs,
                            'Complete verified reference dictionary required',
                        )
                        line = (
                            '\t'.join(
                                columns[:9] + [mapping[s] for s in source_samples],
                            )
                            + '\n'
                        ).encode()
                        ranks = {name: i for i, name in enumerate(declared)}
                        seen_header = True
                else:
                    require(seen_header, 'Missing sample header')
                    fields = text.split('\t')
                    require(len(fields) == 12, 'Exactly three sample columns required')
                    chrom, pos, _, ref, alt = fields[:5]
                    require(chrom in declared and pos.isdigit(), 'Invalid locus')
                    position = int(pos)
                    require(
                        0 < position <= declared[chrom]
                        and position + len(ref) - 1 <= declared[chrom],
                        'Position outside contig',
                    )
                    order = (ranks[chrom], position)
                    require(order >= last, 'Unsorted source')
                    last = order
                    require(
                        re.fullmatch('[ACGTNacgtn]+', ref)
                        and (
                            alt == '.'
                            or all(
                                re.fullmatch('[ACGTNacgtn]+', a) for a in alt.split(',')
                            )
                        ),
                        'Unsupported symbolic/gVCF allele',
                    )
                    require(
                        not any(
                            i.split('=', 1)[0] == 'END' for i in fields[7].split(';')
                        ),
                        'gVCF/END records unsupported',
                    )
                    formats = fields[8].split(':')
                    require(
                        formats[0] == 'GT' and len(formats) == len(set(formats)),
                        'Unique first GT required',
                    )
                    for sample in fields[9:]:
                        values = sample.split(':')
                        require(
                            len(values) <= len(formats)
                            and re.fullmatch(r'(\.|\d+)([/|](\.|\d+))?', values[0]),
                            'Malformed GT/sample',
                        )
                        require(
                            all(
                                int(a) <= (0 if alt == '.' else len(alt.split(',')))
                                for a in re.split('[/|]', values[0])
                                if a != '.'
                            ),
                            'GT outside ALT range',
                        )
                    row_digest.update((text + '\n').encode())
                    records += 1
                dst.write(line)
        require(seen_header, 'Missing header')
        check()
        subprocess.run(  # noqa: S603 - fixed interpreter/code; no shell
            [
                sys.executable,
                '-c',
                'import pysam,sys; pysam.tabix_index(sys.argv[1],preset="vcf",force=False)',
                str(output),
            ],
            check=True,
            timeout=max(0.001, timeout - (time.monotonic() - started)),
            capture_output=True,
        )
        check()
        with pysam.VariantFile(str(output)) as parsed:
            parsed_count = 0
            for record in parsed:
                check()
                require(len(record.samples) == 3, 'Parsed sample count mismatch')
                parsed_count += 1
            require(parsed_count == records, 'Parsed record count mismatch')
        with pysam.TabixFile(str(output)) as indexed:
            count = 0
            index_digest = hashlib.sha256()
            for chrom in indexed.contigs:
                for row in indexed.fetch(chrom):
                    index_digest.update((row + '\n').encode())
                    count += 1
                    check()
            require(
                count == records and index_digest.digest() == row_digest.digest(),
                'Index content/count mismatch',
            )
        require(
            digest(source) == manifest['source_sha256'],
            'Source changed during intake',
        )
        receipt = {
            'status': 'attested_joint_normalized_local_only',
            'prepared': False,
            'searchable': False,
            'clinical_ready': False,
            'manual_load_required': True,
            'family_namespace_integrated': False,
            'calling_history_independently_verified': False,
            'kinship_verified': False,
            'sex_verified': False,
            'manifest': manifest,
            'header_sha256': manifest['header_sha256'],
            'reference_evidence': reference,
            'records': records,
            'sample_order': [mapping[s] for s in source_samples],
            'source_sha256': manifest['source_sha256'],
            'vcf_sha256': digest(output),
            'index_sha256': digest(Path(str(output) + '.tbi')),
            'pysam_version': pysam.__version__,
        }
        (output_dir / 'validation.json').write_text(
            json.dumps(receipt, indent=2) + '\n',
        )
        check()
    except BaseException:
        shutil.rmtree(output_dir)
        raise
    else:
        return receipt


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in (
        'source',
        'output_dir',
        'manifest',
        'reference',
        'family_id',
        'joint_report_id',
    ):
        parser.add_argument(name)
    for name, default in [
        ('max-source-bytes', 536_870_912),
        ('max-expanded-bytes', 536_870_912),
        ('max-line-bytes', 4_194_304),
        ('max-header-bytes', 8_388_608),
        ('max-output-bytes', 536_870_912),
        ('max-scratch-bytes', 1_073_741_824),
        ('timeout', 1800),
    ]:
        parser.add_argument('--' + name, type=int, default=default)
    args = vars(parser.parse_args())
    args['manifest'] = read_json(args['manifest'])
    args['reference'] = read_json(args['reference'])
    sys.stdout.write(json.dumps(normalize(**args), indent=2) + '\n')


if __name__ == '__main__':
    main()
