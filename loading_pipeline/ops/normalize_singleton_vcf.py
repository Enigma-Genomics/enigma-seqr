# ruff: noqa: PLR2004
# Copyright 2026 Enigma Genomics
"""Offline lossless singleton sample rename and fresh BGZF/tabix index (pysam 0.23.3)."""

import argparse
import gzip
import hashlib
import json
import re
import sys
import time
from pathlib import Path

import pysam


def normalize(source, output_dir, expected_sample, sample_id, *, max_source_bytes=2_000_000_000, max_expanded_bytes=536_870_912, max_line_bytes=4_194_304, timeout=1800):  # noqa: C901, PLR0912, PLR0915, PLR0913 - streaming fail-closed checks
    """Write only to a NEW directory; failures leave artifacts for inspection, not use."""
    source, output_dir = Path(source), Path(output_dir)
    if source.stat().st_size > max_source_bytes:
        msg = 'Source exceeds configured byte limit'
        raise ValueError(msg)
    if not isinstance(expected_sample, str) or not re.fullmatch(r'[^\s]{1,200}', expected_sample):
        msg = 'Invalid source sample identifier'
        raise ValueError(msg)
    started = time.monotonic()
    expanded_bytes = 0
    primary_contig_verified = False
    for value in (sample_id,):
        if not isinstance(value, str) or not re.fullmatch(
            r'[A-Za-z0-9][A-Za-z0-9_.-]{0,127}', value,
        ):
            msg = 'Invalid sample identifier'
            raise ValueError(msg)
    output_dir.mkdir(parents=False, exist_ok=False)
    output = output_dir / 'callset.vcf.gz'
    with source.open('rb') as raw:
        compressed = raw.read(2) == b'\x1f\x8b'
    opener = gzip.open if compressed else open
    header_seen, records = False, 0
    reference_seen = False
    with (
        opener(source, 'rt', encoding='utf-8', newline='') as src,
        pysam.BGZFile(str(output), 'w') as dst,
    ):
        while True:
            line = src.readline(max_line_bytes + 1)
            if not line:
                break
            line_size = len(line.encode('utf-8'))
            expanded_bytes += line_size
            if line_size > max_line_bytes or expanded_bytes > max_expanded_bytes or time.monotonic() - started > timeout:
                msg = 'VCF exceeds configured line, expanded-byte or time limit'
                raise ValueError(msg)
            if line.startswith('##contig='):
                match = re.search(r'ID=(?:chr)?1,length=(\d+)(?:,|>)', line)
                if match:
                    if int(match.group(1)) != 248956422:
                        msg = 'Chromosome 1 length is incompatible with GRCh38'
                        raise ValueError(msg)
                    primary_contig_verified = True
            if line.startswith('##reference='):
                reference = line.strip().split('=', 1)[1]
                if reference.rsplit('/', 1)[-1] not in ('GRCh38', 'hg38', 'GRCh38.fa', 'GRCh38.fasta', 'hg38.fa', 'hg38.fasta', 'Homo_sapiens_assembly38.fasta'):
                    msg = 'Requires explicit validated GRCh38 reference'
                    raise ValueError(msg)
                reference_seen = True
            if line.startswith('#CHROM\t'):
                if not reference_seen or not primary_contig_verified:
                    msg = 'Explicit GRCh38 reference and canonical chromosome 1 length required'
                    raise ValueError(msg)
                if header_seen:
                    msg = 'Duplicate VCF header'
                    raise ValueError(msg)
                fields = line.rstrip('\r\n').split('\t')
                if len(fields) != 10 or fields[9] != expected_sample:
                    msg = 'Expected exactly the pinned source sample'
                    raise ValueError(msg)
                fields[9] = sample_id
                line = '\t'.join(fields) + '\n'
                header_seen = True
            elif not line.startswith('#'):
                if not header_seen:
                    msg = 'Missing sample header'
                    raise ValueError(msg)
                fields = line.rstrip('\r\n').split('\t')
                if len(fields) != 10:
                    msg = 'Invalid sample column count'
                    raise ValueError(msg)
                if any(a.startswith('<') or a == '*' or '[' in a or ']' in a for a in fields[4].split(',')):
                    msg = 'Only final SNV/indel VCFs supported; no gVCF or symbolic alleles'
                    raise ValueError(msg)
                formats, values = fields[8].split(':'), fields[9].split(':')
                if formats.count('GT') != 1 or len(set(formats)) != len(formats):
                    msg = 'GT required'
                    raise ValueError(msg)
                gt_index = formats.index('GT')
                if gt_index >= len(values):
                    msg = 'GT missing'
                    raise ValueError(msg)
                gt = values[gt_index]
                if not re.fullmatch(r'(\.|\d+)([/|](\.|\d+))?', gt):
                    msg = 'Unsupported genotype; no rows were silently dropped'
                    raise ValueError(
                        msg,
                    )
                max_allele = 0 if fields[4] == '.' else len(fields[4].split(','))
                if any(int(allele) > max_allele for allele in re.split(r'[/|]', gt) if allele != '.'):
                    msg = 'Genotype allele index outside ALT range'
                    raise ValueError(msg)
                records += 1
            dst.write(line.encode('utf-8'))
    if not header_seen:
        msg = 'Missing sample header'
        raise ValueError(msg)
    pysam.tabix_index(str(output), preset='vcf', force=False)
    with pysam.VariantFile(str(output)) as vcf:
        if list(vcf.header.samples) != [sample_id]:
            msg = 'Output sample mismatch'
            raise ValueError(msg)
        if sum(1 for _ in vcf) != records:
            msg = 'Sequential read count mismatch'
            raise ValueError(msg)
    with pysam.TabixFile(str(output)) as indexed:
        if (
            sum(1 for contig in indexed.contigs for _ in indexed.fetch(contig))
            != records
        ):
            msg = 'Index read count mismatch'
            raise ValueError(msg)

    def digest(path):
        h = hashlib.sha256()
        with path.open('rb') as stream:
            for chunk in iter(lambda: stream.read(1024 * 1024), b''):
                h.update(chunk)
        return h.hexdigest()

    receipt = {
        'status': 'normalized_local_only',
        'prepared': False,
        'searchable': False,
        'source_sha256': digest(source),
        'vcf_sha256': digest(output),
        'index_sha256': digest(Path(str(output) + '.tbi')),
        'records': records,
        'original_sample': expected_sample,
        'sample_id': sample_id,
        'pysam_version': pysam.__version__,
        'manual_load_required': True,
    }
    (output_dir / 'validation.json').write_text(json.dumps(receipt, indent=2) + '\n')
    return receipt


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('source', 'output_dir', 'expected_sample', 'sample_id'):
        parser.add_argument(name)
    args = parser.parse_args()
    sys.stdout.write(json.dumps(normalize(**vars(args)), indent=2) + '\n')


if __name__ == '__main__':
    main()
