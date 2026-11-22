"""Actual ORM adapter tests for the normal isolated Django test database only."""
import hashlib
import json
from unittest.mock import patch

from django.contrib.auth import get_user_model
from django.test import TransactionTestCase

from seqr.models import Individual
from seqr.preparation_models import PreparationLedger, PreparationOutbox
from seqr.utils.preparation_individual import apply_individual_metadata
from seqr.utils.preparation_ledger import apply_preparation, canonical_digest, PreparationConflict


class PreparationIndividualTest(TransactionTestCase):
    databases = {'default'}

    def setUp(self):
        self.user = get_user_model().objects.create(username='synthetic-prep-reviewer', is_superuser=True)
        request = dict(schema_version=1, kind='normalize_singleton', accession_id='SYNTH_PREP', sample_type='WES', genome_version='38',
                       source=dict(bucket='synthetic', path='fixture.vcf', generation='1', vcf_sample_id='OLD'))
        self.envelope = dict(schema_version=2, kind='individual', accession_id=request['accession_id'], sample_type='WES', genome_version='38',
                             source=request['source'], request_id=hashlib.sha256(json.dumps(request, separators=(',', ':')).encode()).hexdigest(),
                             source_sha256='a' * 64, records=1, review_sequence=1, expected_previous_revision=None,
                             status='awaiting_consumer', prepared=False, searchable=False, manual_load_required=True,
                             individual=dict(sample_id='SYNTH_PREP', sex='U', affected='U', metadata=dict(hpo_present=['HP:0000001'], hpo_absent=[])))
        self.seal()
        self.anvil_patch = patch('seqr.models.anvil_enabled', return_value=False)
        self.anvil_patch.start()
        self.addCleanup(self.anvil_patch.stop)

    def seal(self):
        self.envelope['revision'] = canonical_digest({k: v for k, v in self.envelope.items() if k != 'revision'})

    def apply(self):
        return apply_preparation(self.envelope, lambda ledger: apply_individual_metadata(self.envelope, ledger, self.user))

    def next_review(self):
        previous = self.envelope['revision']
        self.envelope.update(review_sequence=2, expected_previous_revision=previous)
        self.envelope['individual']['metadata']['hpo_present'] = []
        self.seal()

    def test_actual_adapter_idempotent_and_explicit_hpo_clear(self):
        first = self.apply()
        self.assertEqual(self.apply().pk, first.pk)
        self.assertEqual(Individual.objects.filter(individual_id='SYNTH_PREP').count(), 1)
        self.next_review()
        self.apply()
        individual = Individual.objects.get(individual_id='SYNTH_PREP')
        self.assertEqual(individual.features, [])
        self.assertEqual(individual.affected, 'U')
        self.assertEqual(PreparationLedger.objects.get(subject_key='SYNTH_PREP').review_sequence, 2)

    def test_manual_edit_blocks_automatic_overwrite(self):
        self.apply()
        individual = Individual.objects.get(individual_id='SYNTH_PREP')
        individual.notes = 'Synthetic manual edit'
        individual.save()
        self.next_review()
        with self.assertRaises(PreparationConflict):
            self.apply()
        individual.refresh_from_db()
        self.assertEqual(individual.notes, 'Synthetic manual edit')
        self.assertEqual(PreparationLedger.objects.get(subject_key='SYNTH_PREP').review_sequence, 1)
        self.assertEqual(PreparationOutbox.objects.count(), 1)

    def test_actual_adapter_rolls_back_new_project_on_failure(self):
        def apply_then_fail(ledger):
            apply_individual_metadata(self.envelope, ledger, self.user)
            raise RuntimeError('Synthetic failure after ORM changes')
        with self.assertRaises(RuntimeError):
            apply_preparation(self.envelope, apply_then_fail)
        self.assertFalse(Individual.objects.filter(individual_id='SYNTH_PREP').exists())
        self.assertFalse(PreparationLedger.objects.filter(subject_key='SYNTH_PREP').exists())
        self.assertFalse(PreparationOutbox.objects.exists())
