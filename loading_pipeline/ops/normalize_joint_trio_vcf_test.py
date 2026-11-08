# Copyright 2026 Enigma Genomics
"""Synthetic-only offline joint intake tests."""

import copy
import gzip
import hashlib
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import pysam

from loading_pipeline.ops.normalize_joint_trio_vcf import (
    manifest_from_request,
    normalize,
)

HEADER = (
    '##fileformat=VCFv4.2\n##reference=GRCh38\n'
    '##contig=<ID=1,length=248956422>\n'
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
    '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tDAD\tCHILD\tMOM\n'
)
ROW = '1\t10\t.\tA\tC\t37.2\tPASS\t.\tGT\t0|1\t./.\t1\n'


class IntakeTest(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.source = self.root / 'input.vcf'
        self.output = self.root / 'out'
        self.reference = {
            'genome_version': '38',
            'verified': True,
            'verification_reference': 'synthetic-reference-fixture',
            'contigs': {'1': 248956422},
        }
        people = []
        for role, source in [
            ('father', 'DAD'),
            ('proband', 'CHILD'),
            ('mother', 'MOM'),
        ]:
            people.append(
                {
                    'role': role,
                    'source_sample_id': source,
                    'sample_id': role,
                    'accession_id': 'acc-' + role,
                    'father_id': 'father' if role == 'proband' else None,
                    'mother_id': 'mother' if role == 'proband' else None,
                },
            )
        self.manifest = {
            'schema_version': 2,
            'kind': 'attested_joint_trio_identity',
            'family_id': 'FAM1',
            'joint_report_id': 'REP1',
            'genome_version': '38',
            'sample_type': 'WGS',
            'dataset_type': 'SNV_INDEL',
            'individuals': people,
            'joint_genotyping_attestation': {
                'confirmed': True,
                'authority': 'user_confirmed',
                'reference': 'synthetic-test-only',
                'confirmed_by': 'test',
            },
        }
        self.write(HEADER + ROW)

    def write(self, text, compressed=False):
        self.source.write_bytes(
            gzip.compress(text.encode()) if compressed else text.encode(),
        )
        header = ''.join(
            line for line in text.splitlines(keepends=True) if line.startswith('#')
        )
        self.manifest['header_sha256'] = hashlib.sha256(header.encode()).hexdigest()
        self.manifest['source_sha256'] = hashlib.sha256(
            self.source.read_bytes(),
        ).hexdigest()

    def run_intake(self, **kwargs):
        return normalize(
            self.source,
            self.output,
            self.manifest,
            self.reference,
            'FAM1',
            'REP1',
            **kwargs,
        )

    def test_lossless_and_fresh_index(self):
        rows = ROW + '1\t11\t.\tA\tC,G\t.\tPASS\t.\tGT\t.|1\t2/2\t.\n'
        self.write(HEADER + rows, compressed=True)
        Path(str(self.source) + '.tbi').write_bytes(b'stale garbage')
        receipt = self.run_intake()
        with gzip.open(self.output / 'callset.vcf.gz', 'rt') as stream:
            self.assertEqual(
                stream.read(),
                HEADER.replace('DAD\tCHILD\tMOM', 'father\tproband\tmother') + rows,
            )
        with pysam.TabixFile(str(self.output / 'callset.vcf.gz')) as indexed:
            self.assertEqual(list(indexed.fetch('1', 10, 11)), [rows.splitlines()[1]])
        self.assertEqual(receipt['records'], 2)
        self.assertFalse(receipt['calling_history_independently_verified'])
        self.assertFalse(receipt['clinical_ready'])

    def test_invalid_records_cleanup(self):
        for row in [
            ROW.replace('0|1', '0/9'),
            ROW.replace('0|1', '0//1'),
            ROW.replace('\tC\t', '\t<NON_REF>\t'),
            ROW.replace('\tPASS\t.\t', '\tPASS\tEND=20\t'),
            ROW.replace('1\t10', '2\t10'),
            ROW.replace('1\t10', '1\t0'),
            ROW + ROW.replace('1\t10', '1\t9'),
            ROW.replace('\t./.\t', '\t./.:5\t'),
        ]:
            with self.subTest(row=row):
                self.write(HEADER + row)
                with self.assertRaises(ValueError):
                    self.run_intake()
                self.assertFalse(self.output.exists())

    def test_invalid_contract(self):
        original = copy.deepcopy(self.manifest)
        for key, value in [
            ('family_id', 'WRONG'),
            ('joint_report_id', 'WRONG'),
            ('genome_version', '37'),
            ('source_sha256', '0' * 64),
            ('joint_genotyping_attestation', {}),
            ('individuals', []),
        ]:
            with self.subTest(key=key):
                self.manifest = dict(original, **{key: value})
                with self.assertRaises(ValueError):
                    self.run_intake()
                self.assertFalse(self.output.exists())
        self.manifest = original
        for field in ('role', 'source_sample_id', 'sample_id', 'accession_id'):
            self.manifest = copy.deepcopy(original)
            self.manifest['individuals'][0][field] = self.manifest['individuals'][1][
                field
            ]
            with self.assertRaises(ValueError):
                self.run_intake()

    def test_header_rejections(self):
        for header in [
            HEADER.replace('MOM', 'OTHER'),
            HEADER.replace('GRCh38', 'GRCh37'),
            HEADER.replace('248956422', '249250621'),
            HEADER + '##reference=GRCh38\n',
            HEADER.replace('##contig=<ID=1,length=248956422>\n', ''),
        ]:
            self.write(header + ROW)
            with self.assertRaises(ValueError):
                self.run_intake()
            self.assertFalse(self.output.exists())

    def test_bounds(self):
        for key in (
            'max_source_bytes',
            'max_expanded_bytes',
            'max_header_bytes',
            'max_line_bytes',
            'max_output_bytes',
            'max_scratch_bytes',
        ):
            with self.subTest(key=key):
                with self.assertRaises(ValueError):
                    self.run_intake(**{key: 10})
                self.assertFalse(self.output.exists())
        with self.assertRaises(ValueError):
            self.run_intake(timeout=0)

    def test_index_failure_and_no_overwrite(self):
        with (
            patch(
                'loading_pipeline.ops.normalize_joint_trio_vcf.subprocess.run',
                side_effect=subprocess.TimeoutExpired('index', 1),
            ),
            self.assertRaises(subprocess.TimeoutExpired),
        ):
            self.run_intake()
        self.assertFalse(self.output.exists())
        self.output.mkdir()
        sentinel = self.output / 'keep'
        sentinel.write_text('preserve')
        with self.assertRaises(FileExistsError):
            self.run_intake()
        self.assertEqual(sentinel.read_text(), 'preserve')

    def test_larger_synthetic_stream(self):
        self.write(
            HEADER + ''.join(ROW.replace('1\t10', f'1\t{i}') for i in range(1, 20001)),
        )
        self.assertEqual(self.run_intake()['records'], 20000)

    def test_identity_only_queue_bridge_and_metadata_rejection(self):
        identity = hashlib.sha256(b'REP1').hexdigest()
        request = {
            'schema_version': 1,
            'kind': 'normalize_joint_trio',
            'source': {'bucket': 'synthetic', 'path': 'input.vcf', 'generation': '1'},
            'subject_key': 'family:' + identity,
            **{
                key: self.manifest[key]
                for key in (
                    'joint_report_id',
                    'genome_version',
                    'sample_type',
                    'dataset_type',
                    'header_sha256',
                    'individuals',
                    'joint_genotyping_attestation',
                )
            },
            'family_id': 'JR_' + identity,
        }

        def seal():
            request['intake_request_id'] = hashlib.sha256(
                json.dumps(
                    {k: v for k, v in request.items() if k != 'intake_request_id'},
                    sort_keys=True,
                    separators=(',', ':'),
                    ensure_ascii=False,
                ).encode(),
            ).hexdigest()

        seal()
        manifest = manifest_from_request(request, self.manifest['source_sha256'])
        receipt = normalize(
            self.source,
            self.output,
            manifest,
            self.reference,
            request['family_id'],
            'REP1',
        )
        self.assertNotIn('metadata', json.dumps(receipt))
        self.assertEqual(receipt['header_sha256'], request['header_sha256'])
        for field, value in [
            ('sex', 'U'),
            ('affected', 'U'),
            ('metadata', {}),
            ('notes', 'test'),
        ]:
            person = request['individuals'][0]
            person[field] = value
            seal()
            with self.assertRaises(ValueError):
                manifest_from_request(request, self.manifest['source_sha256'])
            del person[field]
        request['notes'] = 'not allowed'
        seal()
        with self.assertRaises(ValueError):
            manifest_from_request(request, self.manifest['source_sha256'])

    def test_header_hash_binds_original_crlf_bytes(self):
        self.write(HEADER.replace('\n', '\r\n') + ROW)
        original = self.manifest['header_sha256']
        self.assertNotEqual(original, hashlib.sha256(HEADER.encode()).hexdigest())
        self.assertEqual(self.run_intake()['header_sha256'], original)

    def test_missing_header_newline_rejected(self):
        self.write(HEADER.rstrip('\n'))
        with self.assertRaises(ValueError):
            self.run_intake()
        self.assertFalse(self.output.exists())

    def test_pinned_header_hash_mismatch_cleans_up(self):
        self.manifest['header_sha256'] = '0' * 64
        with self.assertRaises(ValueError):
            self.run_intake()
        self.assertFalse(self.output.exists())

    def test_cli(self):
        manifest = self.root / 'manifest.json'
        reference = self.root / 'reference.json'
        manifest.write_text(json.dumps(self.manifest))
        reference.write_text(json.dumps(self.reference))
        script = Path(__file__).with_name('normalize_joint_trio_vcf.py')
        result = subprocess.run(  # noqa: S603 - synthetic fixed CLI
            [
                sys.executable,
                str(script),
                str(self.source),
                str(self.output),
                str(manifest),
                str(reference),
                'FAM1',
                'REP1',
                '--timeout',
                '60',
            ],
            check=True,
            capture_output=True,
            text=True,
        )
        self.assertEqual(json.loads(result.stdout)['records'], 1)


if __name__ == '__main__':
    unittest.main()
