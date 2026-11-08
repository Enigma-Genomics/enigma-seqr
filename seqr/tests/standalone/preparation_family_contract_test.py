# ruff: noqa: E402 - standalone Django configuration precedes model imports

"""Standalone synthetic contract tests; no application or clinical database."""

import copy
import hashlib
import unittest

from django.conf import settings

settings.configure(
    INSTALLED_APPS=[],
    DATABASES={'default': {'ENGINE': 'django.db.backends.sqlite3', 'NAME': ':memory:'}},
    USE_TZ=True,
)
import django

django.setup()

from seqr.utils.preparation_family_contract import (
    family_id,
    family_subject_key,
    family_request,
    validate_family_envelope,
)
from seqr.utils.preparation_family_test_utils import seal, synthetic_family_envelope
from seqr.utils.preparation_ledger import canonical_digest, PreparationConflict


class FamilyContractTest(unittest.TestCase):
    def setUp(self):
        self.envelope = synthetic_family_envelope()

    def test_valid_explicit_order_and_roles(self):
        self.assertEqual(
            validate_family_envelope(self.envelope), self.envelope['subject_key']
        )
        request = family_request(self.envelope)
        self.envelope['individuals'].reverse()
        seal(self.envelope)
        self.assertEqual(family_request(self.envelope), request)
        validate_family_envelope(self.envelope)

    def test_maximum_identity_cross_language_rule(self):
        value = 'X' * 128
        expected = hashlib.sha256(value.encode('utf-8')).hexdigest()
        self.assertEqual(family_id(value), 'JR_' + expected)
        self.assertEqual(family_subject_key(value), 'family:' + expected)
        self.assertEqual(len(family_id(value)), 67)
        self.assertEqual(len(family_subject_key(value)), 71)
        for bad in ('', 'X' * 129, 'family:123', ' leading', '日本語'):
            with self.assertRaises(PreparationConflict):
                family_id(bad)

    def test_metadata_revision_preserves_request(self):
        prior = self.envelope['request_id']
        self.envelope['individuals'][0]['sex'] = 'F'
        self.envelope['individuals'][0]['metadata']['hpo_present'] = []
        self.envelope['metadata_review']['reference'] = 'synthetic-review-2'
        seal(self.envelope)
        validate_family_envelope(self.envelope)
        self.assertEqual(canonical_digest(family_request(self.envelope)), prior)

    def test_immutable_request_covers_source_member_role_assay(self):
        for mutate in (
            lambda e: e['source'].update(generation='2'),
            lambda e: e.update(source_sha256='d' * 64),
            lambda e: e.update(sample_type='WGS'),
            lambda e: e['individuals'][0].update(accession_id='OTHER'),
            lambda e: e['individuals'][0].update(source_sample_id='OTHER'),
            lambda e: e['joint_genotyping_attestation'].update(reference='other'),
        ):
            envelope = copy.deepcopy(self.envelope)
            mutate(envelope)
            seal(envelope)
            with self.assertRaises(PreparationConflict):
                validate_family_envelope(envelope)

    def test_invalid_metadata_pedigree_and_claims(self):
        for mutate in (
            lambda e: e.update(family_id='JR_wrong'),
            lambda e: e.update(subject_key='SYNTHETIC_JOINT_REPORT'),
            lambda e: e.update(clinical_ready=True),
            lambda e: e.update(calling_history_independently_verified=True),
            lambda e: e['metadata_review'].update(confirmed=False),
            lambda e: e['individuals'][0].update(role='mother'),
            lambda e: e['individuals'][0].update(mother_id='SYNTH_mother'),
            lambda e: e['individuals'][0].update(sex='inferred'),
            lambda e: e['individuals'][0]['metadata'].update(hpo_absent=['HP:0000001']),
            lambda e: e['individuals'][0]['metadata'].update(notes='not yet supported'),
        ):
            envelope = copy.deepcopy(self.envelope)
            mutate(envelope)
            seal(envelope)
            with self.assertRaises(PreparationConflict):
                validate_family_envelope(envelope)


if __name__ == '__main__':
    unittest.main()
