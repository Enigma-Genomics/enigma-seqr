# Copyright 2026 Enigma contributors
"""Bounded synthetic-only local resource smoke test; not a clinical/cloud benchmark."""

import argparse
import hashlib
import json
import resource
import shutil
import tempfile
import time
from pathlib import Path

import pysam

from loading_pipeline.ops.merge_trio_vcfs import ROLES, merge

HEADER = """##fileformat=VCFv4.2
##reference=GRCh38
##contig=<ID=chr1,length=10000000>
##INFO=<ID=TEST_NOTE,Number=1,Type=String,Description="Synthetic incompressible workload">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allele depth">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Likelihood">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSYNTHETIC
"""


def benchmark(rows):
    if not 1 <= rows <= 1000000:  # noqa: PLR2004 - explicit synthetic safety bound
        msg = 'Choose 1 to 1000000 synthetic rows per source'
        raise ValueError(msg)
    with tempfile.TemporaryDirectory(prefix='synthetic-trio-benchmark-') as tmp:
        root = Path(tmp)
        first = root / 'proband.vcf.gz'
        with pysam.BGZFile(str(first), 'w') as stream:
            stream.write(HEADER.encode())
            for pos in range(1, rows + 1):
                note = hashlib.shake_256(str(pos).encode()).hexdigest(256)
                stream.write(
                    f'chr1\t{pos}\t.\tA\tC\t.\tPASS\tTEST_NOTE={note}\tGT:AD:PL\t0/1:8,4:10,0,20\n'.encode(),
                )
        members = []
        for role in ROLES:
            path = root / f'{role}.vcf.gz'
            if path != first:
                shutil.copyfile(first, path)
            members.append(
                {
                    'role': role,
                    'path': str(path),
                    'sample_id': role,
                    'vcf_sample_id': 'SYNTHETIC',
                    'genome_version': 'GRCh38',
                },
            )
        start = time.monotonic()
        receipt = merge(
            {'genome_version': 'GRCh38', 'members': members},
            root / 'output',
        )
        return {
            'synthetic_only': True,
            'rows_per_source': rows,
            'compressed_input_bytes': 3 * first.stat().st_size,
            'expanded_input_bytes': receipt['expanded_source_bytes'],
            'output_records': receipt['records'],
            'scratch_peak_sampled_bytes': receipt['scratch_peak_sampled_bytes'],
            'wall_seconds': round(time.monotonic() - start, 3),
            'child_maxrss_platform_units': resource.getrusage(
                resource.RUSAGE_CHILDREN,
            ).ru_maxrss,
            'rss_units': 'bytes on macOS; KiB on Linux',
            'scope': 'All loci shared; identical scalar annotation; not a GATK compatibility or cloud-cost benchmark',
        }


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--rows', type=int, default=10000)
    print(json.dumps(benchmark(parser.parse_args().rows), indent=2))  # noqa: T201 - explicit benchmark CLI output
