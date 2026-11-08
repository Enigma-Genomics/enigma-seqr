# Copyright 2026 Enigma Genomics
import gzip
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from loading_pipeline.ops.normalize_singleton_vcf import normalize

HEADER = (
    '##fileformat=VCFv4.2\n##reference=GRCh38\n##contig=<ID=chr1,length=1000>\n'
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
    '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">\n'
    '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\told\n'
)
BODY = (
    'chr1\t10\t.\tA\tT\t.\tPASS\t.\tGT:DP\t0|1:20\n'
    'chr1\t20\t.\tG\tC\t.\tPASS\t.\tGT:DP\t./.:.\n'
    'chr1\t30\t.\tC\tA,G\t.\tPASS\t.\tGT:DP\t1/2:30\n'
)


class NormalizeTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.src = self.root / 'source.vcf'
        self.src.write_text(HEADER + BODY)
        self.out = self.root / 'output'

    def test_lossless_rename_and_fresh_index(self):
        receipt = normalize(self.src, self.out, 'old', 'new')
        self.assertEqual(receipt['records'], 3)
        with gzip.open(self.out / 'callset.vcf.gz', 'rt') as stream:
            self.assertEqual(stream.read(), HEADER.replace('\told\n', '\tnew\n') + BODY)
        self.assertTrue((self.out / 'callset.vcf.gz.tbi').exists())
        self.assertFalse(receipt['searchable'])
        with self.assertRaises(FileExistsError):
            normalize(self.src, self.out, 'old', 'new')

    def test_wrong_sample_no_receipt(self):
        with self.assertRaises(ValueError):
            normalize(self.src, self.out, 'wrong', 'new')
        self.assertFalse((self.out / 'validation.json').exists())

    def test_polyploid_not_dropped(self):
        self.src.write_text(HEADER + BODY.replace('1/2:30', '0/1/2:30'))
        with self.assertRaises(ValueError):
            normalize(self.src, self.out, 'old', 'new')
        self.assertFalse((self.out / 'validation.json').exists())

    def test_unsorted_no_receipt(self):
        self.src.write_text(
            HEADER + '\n'.join(reversed(BODY.strip().splitlines())) + '\n',
        )
        with self.assertRaises(OSError):
            normalize(self.src, self.out, 'old', 'new')
        self.assertFalse((self.out / 'validation.json').exists())

    def test_reference_and_symbolic_fail_closed(self):
        for index, content in enumerate([
            (HEADER + BODY).replace('##reference=GRCh38\n', ''),
            (HEADER + BODY).replace('GRCh38', 'GRCh37'),
            HEADER + BODY.replace('A,G', '<NON_REF>'),
        ]):
            self.src.write_text(content)
            with self.subTest(index=index), self.assertRaises(ValueError):
                normalize(self.src, self.root / f'out-{index}', 'old', 'new')

    def test_invalid_duplicate_gt_and_allele_index(self):
        for index, content in enumerate([
            HEADER + BODY.replace('GT:DP', 'GT:GT'),
            HEADER + BODY.replace('0|1:20', '0/2:20'),
        ]):
            self.src.write_text(content)
            with self.subTest(index=index), self.assertRaises(ValueError):
                normalize(self.src, self.root / f'invalid-{index}', 'old', 'new')

    def test_cli_executes(self):
        result = subprocess.run(  # noqa: S603 - trusted local module and synthetic temp paths
            [sys.executable, '-m', 'loading_pipeline.ops.normalize_singleton_vcf',
             str(self.src), str(self.out), 'old', 'new'],
            check=True, capture_output=True, text=True,
        )
        self.assertEqual(json.loads(result.stdout)['records'], 3)
        self.assertTrue((self.out / 'validation.json').exists())
