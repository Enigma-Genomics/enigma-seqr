"""Actual ORM adapter tests for the normal isolated Django test database only."""
import copy
import hashlib
import json
from unittest.mock import patch

from django.contrib.auth import get_user_model
from django.test import TransactionTestCase
from django.utils import timezone

from seqr.models import Individual, Project, Dataset
from seqr.preparation_models import PreparationLedger, PreparationOutbox
from seqr.utils.preparation_individual import apply_individual_metadata, refresh_preparation_dictionaries
from seqr.utils.preparation_ledger import apply_preparation, canonical_digest, PreparationConflict, validate_envelope, deliver_acknowledgment


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

    def test_reject_malformed_envelopes_before_database_changes(self):
        original = copy.deepcopy(self.envelope)
        mutations = [
            ('schema', lambda e: e.update(schema_version=99)),
            ('identity', lambda e: e.update(accession_id='bad/id')),
            ('reserved', lambda e: e.update(accession_id='JR_' + 'a' * 64)),
            ('request', lambda e: e.update(request_id='invalid')),
            ('source', lambda e: e['source'].update(generation='0')),
            ('fingerprint', lambda e: e.update(request_id='b' * 64)),
            ('normalization', lambda e: e.update(records=True)),
            ('sequence', lambda e: e.update(review_sequence=True)),
            ('predecessor', lambda e: e.update(expected_previous_revision='a' * 64)),
            ('person', lambda e: e['individual'].update(sex='invalid')),
            ('metadata', lambda e: e['individual']['metadata'].update(unknown=True)),
            ('hpo', lambda e: e['individual']['metadata'].update(hpo_present=['invalid'])),
            ('hpo conflict', lambda e: e['individual']['metadata'].update(hpo_absent=['HP:0000001'])),
            ('year', lambda e: e['individual']['metadata'].update(birth_year=True)),
            ('boolean', lambda e: e['individual']['metadata'].update(consanguinity='false')),
            ('text', lambda e: e['individual']['metadata'].update(notes=1)),
            ('onset', lambda e: e['individual']['metadata'].update(onset_age='invalid')),
            ('inheritance', lambda e: e['individual']['metadata'].update(expected_inheritance=['invalid'])),
            ('lifecycle', lambda e: e.update(searchable=True)),
        ]
        for name, mutate in mutations:
            with self.subTest(name=name):
                self.envelope = copy.deepcopy(original)
                mutate(self.envelope)
                self.seal()
                with self.assertRaises(PreparationConflict):
                    validate_envelope(self.envelope)
        self.envelope = copy.deepcopy(original)
        self.envelope['revision'] = '0' * 64
        with self.assertRaises(PreparationConflict):
            validate_envelope(self.envelope)
        self.assertFalse(PreparationLedger.objects.exists())

    def test_reject_self_consistent_unsupported_build(self):
        self.envelope['genome_version'] = '37'
        request = dict(schema_version=1, kind='normalize_singleton', **{
            k: self.envelope[k] for k in ('accession_id', 'sample_type', 'genome_version', 'source')})
        self.envelope['request_id'] = hashlib.sha256(json.dumps(request, separators=(',', ':')).encode()).hexdigest()
        self.seal()
        with self.assertRaisesRegex(PreparationConflict, 'dataset provenance'):
            self.apply()

    def test_ordered_cas_invalid_adapter_and_ack_recovery(self):
        with self.assertRaisesRegex(PreparationConflict, 'state token'):
            apply_preparation(self.envelope, lambda ledger: ('R_INVALID', 'bad'))
        self.assertFalse(PreparationLedger.objects.exists())
        first = self.apply()
        self.next_review()
        self.envelope['review_sequence'] = 3
        self.seal()
        with self.assertRaisesRegex(PreparationConflict, 'out-of-order'):
            self.apply()
        sent = []
        def fail_after_publish(payload):
            sent.append(payload)
            raise RuntimeError('Synthetic post-upload failure')
        with self.assertRaises(RuntimeError):
            deliver_acknowledgment(first.pk, fail_after_publish)
        first.refresh_from_db()
        self.assertIsNone(first.delivered_at)
        deliver_acknowledgment(first.pk, sent.append)
        deliver_acknowledgment(first.pk, sent.append)
        self.assertEqual(sent, [first.payload, first.payload])

    def test_supported_metadata_mapping_and_dictionary_refresh(self):
        self.envelope['individual']['metadata'].update(birth_year=2001, death_year=2020,
            notes='Synthetic note', consanguinity=False, affected_relatives=True,
            onset_age='', expected_inheritance=['R'], candidate_genes='')
        self.seal()
        self.apply()
        individual = Individual.objects.get()
        self.assertEqual(individual.birth_year, 2001)
        self.assertEqual(individual.death_year, 2020)
        self.assertEqual(individual.notes, 'Synthetic note')
        self.assertIsNone(individual.onset_age)
        self.assertEqual(individual.expected_inheritance, ['R'])
        self.assertEqual(individual.candidate_genes, [])
        with patch('clickhouse_search.models.postgres_dicts.AffectedDict.reload') as affected, patch('clickhouse_search.models.postgres_dicts.SexDict.reload') as sex, patch('clickhouse_search.models.postgres_dicts.IndividualMetadataDict.reload') as metadata:
            refresh_preparation_dictionaries(self.user)
            for mock in (affected, sex, metadata):
                mock.assert_called_once_with(self.user)

    def test_adapter_direct_identity_and_free_text_rejection(self):
        ledger = PreparationLedger(subject_key='SYNTH_PREP')
        self.envelope['individual']['sample_id'] = 'OTHER'
        with self.assertRaisesRegex(PreparationConflict, 'identity or reference'):
            apply_individual_metadata(self.envelope, ledger, self.user)
        self.envelope['individual']['sample_id'] = 'SYNTH_PREP'
        self.envelope['individual']['metadata']['candidate_genes'] = 'unmapped free text'
        self.seal()
        with self.assertRaisesRegex(PreparationConflict, 'structured mapping'):
            self.apply()

    def test_reject_ambiguous_project(self):
        self.apply()
        Project.objects.create(name='SYNTH_PREP', genome_version='38', created_by=self.user)
        self.next_review()
        with self.assertRaisesRegex(PreparationConflict, 'Ambiguous'):
            self.apply()

    def test_reject_project_reference_and_identity_change(self):
        self.apply()
        project = Project.objects.get()
        project.genome_version = '37'
        project.save()
        self.next_review()
        with self.assertRaisesRegex(PreparationConflict, 'Project reference'):
            self.apply()

    def test_reject_changed_singleton_membership(self):
        self.apply()
        Individual.objects.update(individual_id='OTHER')
        self.next_review()
        with self.assertRaisesRegex(PreparationConflict, 'expected singleton'):
            self.apply()

    def test_reject_unreviewed_existing_record_adoption(self):
        self.apply()
        PreparationOutbox.objects.all().delete()
        PreparationLedger.objects.all().delete()
        with self.assertRaisesRegex(PreparationConflict, 'baseline adoption'):
            self.apply()

    def test_reject_loaded_and_missing_project(self):
        self.apply()
        dataset = Dataset.objects.create(sample_type='WES', dataset_type=Dataset.DATASET_TYPE_VARIANT_CALLS,
            data_source='synthetic', loaded_date=timezone.now(), created_by=self.user)
        individual = Individual.objects.get()
        self.next_review()
        for relation in (dataset.active_individuals, dataset.inactive_individuals):
            relation.add(individual)
            with self.assertRaisesRegex(PreparationConflict, 'Loaded individual'):
                self.apply()
            relation.clear()
        Project.objects.get().delete()
        with self.assertRaisesRegex(PreparationConflict, 'project is missing'):
            self.apply()
