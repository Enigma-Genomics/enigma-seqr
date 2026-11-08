# Copyright 2026 Enigma contributors
# ruff: noqa: C901, PLR0912, PLR0915, PLR2004, PLW2901
"""Fail-closed offline final-VCF union. This is NOT joint genotyping or clinical QC."""

import argparse
import gzip
import hashlib
import json
import math
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import pysam
import pysam.bcftools

ROLES = ('proband', 'mother', 'father')


def require(condition, message):
    if not condition:
        raise ValueError(message)


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda: stream.read(1048576), b''):
            h.update(chunk)
    return h.hexdigest()


def vectors(alleles, ploidy):
    if ploidy == 1:
        return [(a,) for a in alleles]
    return [(alleles[i], alleles[j]) for j in range(len(alleles)) for i in range(j + 1)]


def check_call(source, result):
    """Verify each source FORMAT value after allele remapping, including partial GT."""
    original = source.samples[0]
    merged = result.samples[original.name]
    require(original.phased == merged.phased, 'Phase changed')
    for key, value in original.items():
        require(key in merged, f'FORMAT dropped: {key}')
        got = merged[key]
        number = source.header.formats[key].number
        if key == 'GT':
            wanted = tuple(None if a is None else source.alleles[a] for a in value)
            actual = tuple(None if a is None else result.alleles[a] for a in got)
            require(wanted == actual, 'GT alleles changed')
        elif number in ('A', 'R', 'G') and value != (None,):
            old = source.alleles[1:] if number == 'A' else source.alleles
            new = result.alleles[1:] if number == 'A' else result.alleles
            if number == 'G':
                old, new = (
                    vectors(old, len(original['GT'])),
                    vectors(new, len(original['GT'])),
                )
                old, new = [frozenset(x) for x in old], [frozenset(x) for x in new]
            require(
                all(value[i] == got[new.index(a)] for i, a in enumerate(old)),
                f'FORMAT remap changed: {key}',
            )
        else:
            require(value == got, f'FORMAT changed: {key}')


def _build(manifest, staging, max_records, max_bytes):
    require(pysam.__version__ == '0.23.3', 'Vetted pysam 0.23.3 required')
    require(
        manifest.get('genome_version') in ('GRCh37', 'GRCh38'),
        'Explicit genome_version required',
    )
    members = manifest.get('members', [])
    require(
        len(members) == 3 and {m['role'] for m in members} == set(ROLES),
        'Exactly three explicit unique roles required',
    )
    members = sorted(members, key=lambda m: ROLES.index(m['role']))
    targets = [m['sample_id'] for m in members]
    require(len(set(targets)) == 3, 'Duplicate target IDs')
    for m in members:
        require(
            all(
                re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]{0,127}', m[k])
                for k in ('sample_id', 'vcf_sample_id')
            ),
            'Invalid sample ID',
        )
        require(
            m.get('genome_version') == manifest['genome_version'],
            'Source genome mismatch',
        )
    paths = [Path(m['path']).resolve(strict=True) for m in members]
    require(len(set(paths)) == 3, 'Duplicate source paths')
    require(
        sum(p.stat().st_size for p in paths) <= max_bytes,
        'Input byte limit exceeded',
    )
    hashes = [digest(p) for p in paths]
    for member, sha in zip(members, hashes, strict=True):
        require(
            member.get('expected_sha256', sha) == sha, 'Pinned source SHA256 mismatch',
        )
    sources, records, contigs, definitions, site_metadata = [], {}, None, {}, {}
    counts = []
    expanded_bytes = 0
    for n, (m, path) in enumerate(zip(members, paths, strict=True)):
        renamed = staging / f'{n}.vcf.gz'
        with path.open('rb') as stream:
            opener = gzip.open if stream.read(2) == b'\x1f\x8b' else open
        # Preserve all body/header fields, rename only the explicitly pinned sample.
        with opener(path, 'rt') as src, pysam.BGZFile(str(renamed), 'w') as dst:
            seen, line_count = False, 0
            for line in src:
                expanded_bytes += len(line.encode())
                require(
                    expanded_bytes <= 536870912,
                    'Expanded input limit exceeded (512 MiB)',
                )
                require(len(line) <= 1048576, 'VCF line limit exceeded')
                if line.startswith('#CHROM\t'):
                    f = line.rstrip('\r\n').split('\t')
                    require(
                        not seen and len(f) == 10 and f[9] == m['vcf_sample_id'],
                        'Exact singleton source sample mismatch',
                    )
                    f[9], seen = m['sample_id'], True
                    line = '\t'.join(f) + '\n'
                elif not line.startswith('#'):
                    require(seen, 'Missing sample header')
                    f = line.rstrip('\r\n').split('\t')
                    require(len(f) == 10, 'Invalid sample columns')
                    keys, values = f[8].split(':'), f[9].split(':')
                    require(
                        len(keys) == len(set(keys))
                        and len(keys) == len(values)
                        and 'GT' in keys,
                        'Invalid FORMAT',
                    )
                    require(
                        re.fullmatch(
                            r'(\.|\d+)([/|](\.|\d+))?',
                            values[keys.index('GT')],
                        ),
                        'Unsupported GT',
                    )
                    gt_index = keys.index('GT')
                    order = [gt_index] + [i for i in range(len(keys)) if i != gt_index]
                    f[8] = ':'.join(keys[i] for i in order)
                    f[9] = ':'.join(values[i] for i in order)
                    line = '\t'.join(f) + '\n'
                    line_count += 1
                    require(line_count <= max_records, 'Record limit exceeded')
                dst.write(line.encode())
            require(seen, 'Missing sample header')
        with pysam.VariantFile(str(renamed)) as vcf:
            refs = [r.value for r in vcf.header.records if r.key == 'reference']
            require(
                refs == [manifest['genome_version']],
                'Exact ##reference=genome_version required; no build inference',
            )
            current = tuple((k, v.length) for k, v in vcf.header.contigs.items())
            require(
                current and all(v is not None and v > 0 for _, v in current),
                'Explicit contig lengths required',
            )
            require(contigs is None or current == contigs, 'Contig dictionary mismatch')
            contigs = current
            for category in ('formats', 'info', 'filters'):
                for key, entry in getattr(vcf.header, category).items():
                    definition = str(entry.record)
                    require(
                        (category, key) not in definitions
                        or definitions[category, key] == definition,
                        'Conflicting header definition',
                    )
                    definitions[category, key] = definition
            last, count = (-1, -1), 0
            for row in vcf:
                order = (list(vcf.header.contigs).index(row.contig), row.pos)
                require(
                    order > last,
                    'Unsorted or duplicate source locus; normalize separately',
                )
                last = order
                require(
                    row.stop == row.start + len(row.ref)
                    and 'END=' not in str(row).split('\t')[7],
                    'gVCF/reference blocks unsupported: joint genotype separately',
                )
                require(
                    row.alts
                    and all(re.fullmatch('[ACGTNacgtn]+', a) for a in row.alleles),
                    'Symbolic/gVCF/breakend alleles unsupported',
                )
                require(
                    len(row.alleles) <= 32
                    and row.stop <= vcf.header.contigs[row.contig].length,
                    'Allele/contig bounds exceeded',
                )
                call = row.samples[0]
                gt = call['GT']
                require(
                    all(a is None or 0 <= a < len(row.alleles) for a in gt),
                    'GT allele out of range',
                )
                for key, value in call.items():
                    require(('formats', key) in definitions, 'Undefined FORMAT')
                    number = vcf.header.formats[key].number
                    if key != 'GT' and value != (None,) and number in ('A', 'R', 'G'):
                        length = {
                            'A': len(row.alleles) - 1,
                            'R': len(row.alleles),
                            'G': math.comb(len(row.alleles) + len(gt) - 1, len(gt)),
                        }[number]
                        require(
                            len(value) == length,
                            'Malformed allele-vector cardinality',
                        )
                locus = (row.contig, row.pos)
                # Scalar INFO must agree at shared loci; allele-dependent INFO cannot be safely combined here.
                require(
                    all(
                        vcf.header.info[k].number not in ('A', 'R', 'G')
                        for k in row.info
                    ),
                    'Allele-dependent INFO unsupported; do not drop it, use vetted joint workflow',
                )
                require(
                    all(('info', k) in definitions for k in row.info),
                    'Undefined INFO',
                )
                require(
                    all(('filters', k) in definitions for k in row.filter),
                    'Undefined FILTER',
                )
                metadata = (
                    row.ref,
                    row.id,
                    row.qual,
                    tuple(row.filter),
                    dict(row.info),
                )
                require(
                    locus not in site_metadata or site_metadata[locus] == metadata,
                    'Conflicting shared-site INFO/REF/ID/QUAL/FILTER; no lossy reconciliation',
                )
                site_metadata[locus] = metadata
                records[m['sample_id'], locus] = row.copy()
                count += 1
            counts.append(count)
        pysam.tabix_index(str(renamed), preset='vcf')  # Never trust an input index.
        sources.append(str(renamed))
    output = staging / 'callset.vcf.gz'
    pysam.bcftools.merge(
        '--no-version',
        '-m',
        'all',
        '-i',
        '-',
        '-O',
        'z',
        '-o',
        str(output),
        *sources,
        catch_stdout=False,
    )
    pysam.tabix_index(str(output), preset='vcf')
    visited = set()
    with pysam.VariantFile(str(output)) as vcf:
        require(list(vcf.header.samples) == targets, 'Merged sample mismatch')
        for row in vcf:
            locus = (row.contig, row.pos)
            require(
                locus not in visited and locus in site_metadata,
                'Output duplicated or added locus',
            )
            visited.add(locus)
            require(
                (row.ref, row.id, row.qual, tuple(row.filter), dict(row.info))
                == site_metadata[locus],
                'Site metadata changed',
            )
            for sample in targets:
                source = records.get((sample, locus))
                if source is None:
                    require(
                        all(a is None for a in row.samples[sample]['GT']),
                        'Absent site became a called genotype',
                    )
                else:
                    check_call(source, row)
    require(visited == set(site_metadata), 'Source locus dropped')
    with pysam.TabixFile(str(output)) as indexed:
        require(
            sum(1 for c in indexed.contigs for _ in indexed.fetch(c)) == len(visited),
            'Fresh index mismatch',
        )
    require(hashes == [digest(p) for p in paths], 'Source changed during preparation')
    receipt = {
        'status': 'merged_local_only',
        'joint_genotyped': False,
        'prepared': False,
        'searchable': False,
        'manual_load_required': True,
        'genome_version': manifest['genome_version'],
        'pysam_version': pysam.__version__,
        'members': [
            {**m, 'sha256': h, 'records': c}
            for m, h, c in zip(members, hashes, counts, strict=True)
        ],
        'records': len(visited),
        'vcf_sha256': digest(output),
        'index_sha256': digest(str(output) + '.tbi'),
        'qc_limits': [
            'No kinship, sex or clinical validation',
            'Not joint genotyped',
            'Exomiser trio unsupported until separately tested',
        ],
    }
    (staging / 'validation.json').write_text(json.dumps(receipt, indent=2) + '\n')
    for source in sources:
        Path(source).unlink()
        Path(source + '.tbi').unlink()
    return receipt


def merge(
    manifest,
    output_dir,
    timeout=1800,
    max_records=1000000,
    max_bytes=2000000000,
):
    """Subprocess bounds wall time; byte/record/line/allele limits bound work. No network."""
    output = Path(output_dir).absolute()
    require(not output.exists(), 'Output must not exist')
    require(
        timeout > 0 and 0 < max_records <= 1000000 and 0 < max_bytes <= 2000000000,
        'Invalid resource limits',
    )
    # Reserve output exclusively, never overwrite or touch an existing destination.
    output.mkdir()
    try:
        with tempfile.TemporaryDirectory(prefix='.trio-', dir=output.parent) as tmp:
            staging = Path(tmp)
            request = staging / 'request.json'
            request.write_text(json.dumps(manifest))
            result = subprocess.run(  # noqa: S603 - fixed local worker, no shell
                [
                    sys.executable,
                    str(Path(__file__).resolve()),
                    '--worker',
                    str(request),
                    str(staging),
                    str(max_records),
                    str(max_bytes),
                ],
                capture_output=True,
                text=True,
                timeout=timeout,
                check=False,
            )
            require(
                result.returncode == 0,
                'Trio preparation failed: ' + result.stderr[-4000:],
            )
            for name in ('callset.vcf.gz', 'callset.vcf.gz.tbi', 'validation.json'):
                shutil.move(str(staging / name), str(output / name))
        return json.loads((output / 'validation.json').read_text())
    except BaseException:
        shutil.rmtree(output)
        raise


def main():
    if len(sys.argv) > 1 and sys.argv[1] == '--worker':
        _build(
            json.loads(Path(sys.argv[2]).read_text()),
            Path(sys.argv[3]),
            int(sys.argv[4]),
            int(sys.argv[5]),
        )
        return
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('manifest')
    parser.add_argument('output_dir')
    parser.add_argument('--timeout', type=int, default=1800)
    args = parser.parse_args()
    merge(
        json.loads(Path(args.manifest).read_text()),
        args.output_dir,
        timeout=args.timeout,
    )


if __name__ == '__main__':
    main()
