# Copyright 2026 Enigma Genomics
import copy
import unittest

from loading_pipeline.ops.family_prep_plan import ManifestError, prep_plan


class FamilyPlanTest(unittest.TestCase):
    def setUp(self):
        people, observations = [], {}
        for role, sid in [('proband', 'CHILD'), ('father', 'DAD'), ('mother', 'MOM')]:
            source = {
                'bucket': 'synthetic',
                'path': f'{sid}.vcf.gz',
                'generation': '1',
                'vcf_sample_id': sid,
                'genome_version': '38',
                'sample_type': 'WES',
                'dataset_type': 'SNV_INDEL',
            }
            people.append(
                {
                    'role': role,
                    'accession_id': sid,
                    'sample_id': sid,
                    'source': source,
                    'father_id': 'DAD' if role == 'proband' else None,
                    'mother_id': 'MOM' if role == 'proband' else None,
                    'sex': 'U',
                    'affected': 'U',
                    'metadata': {},
                },
            )
            observations[sid] = dict(source, sample_ids=[sid])
        people[0]['metadata'] = {
            'hpo_present': ['HP:0000001'], 'onset_age': 'C', 'expected_inheritance': ['R'],
        }
        self.m = {
            'schema_version': 2,
            'kind': 'family',
            'family_id': 'SYNTHFAMILY',
            'joint_report_id': 'SYNTHREPORT',
            'genome_version': '38',
            'sample_type': 'WES',
            'manual_load_required': True,
            'individuals': people,
        }
        self.obs = observations

    def test_per_person_mapping(self):
        result = prep_plan(self.m, self.obs)
        self.assertFalse(result['prepared'])
        self.assertEqual(result['individuals'][0]['father_id'], 'DAD')
        self.assertEqual(result['individuals'][0]['onset_age'], 'C')
        self.assertEqual(result['individuals'][1]['features'], [])
        self.assertEqual(result['individuals'][1]['affected'], 'U')

    def test_fail_closed_sources(self):
        for field, value in [
            ('generation', '2'),
            ('sample_ids', ['OTHER']),
            ('genome_version', '37'),
            ('sample_type', 'WGS'),
            ('dataset_type', 'SV'),
        ]:
            with self.subTest(field=field), self.assertRaises(ManifestError):
                changed = copy.deepcopy(self.obs)
                changed['CHILD'][field] = value
                prep_plan(self.m, changed)

    def test_invalid_pedigree_and_phenotypes(self):
        for field, value in [
            ('father_id', 'MOM'),
            ('sample_id', 'DAD'),
            ('affected', 'unknown'),
            ('metadata', {'hpo_present': ['HP:0000001'], 'hpo_absent': ['HP:0000001']}),
        ]:
            with self.subTest(field=field), self.assertRaises(ManifestError):
                changed = copy.deepcopy(self.m)
                changed['individuals'][0][field] = value
                prep_plan(changed, self.obs)
