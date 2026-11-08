"""Actual family ORM tests. Run only using Django's disposable test database."""

import copy
from unittest.mock import patch

from django.contrib.auth import get_user_model
from django.test import TransactionTestCase, override_settings
from django.utils import timezone

from seqr.models import Project, Family, Individual, Dataset
from seqr.preparation_models import PreparationLedger, PreparationOutbox
from seqr.utils.preparation_family import apply_family_metadata
from seqr.utils.preparation_family_contract import family_request
from seqr.utils.preparation_family_test_utils import synthetic_family_envelope, seal
from seqr.utils.preparation_ledger import (
    apply_preparation,
    PreparationConflict,
    canonical_digest,
)


@override_settings(PREPARATION_FAMILY_ENABLED=True)
class PreparationFamilyTest(TransactionTestCase):
    databases = {'default'}

    def setUp(self):
        self.user = get_user_model().objects.create(
            username='synthetic-family-reviewer', is_superuser=True
        )
        self.envelope = synthetic_family_envelope()
        self.anvil = patch('seqr.models.anvil_enabled', return_value=False)
        self.anvil.start()
        self.addCleanup(self.anvil.stop)

    def apply(self):
        return apply_preparation(
            self.envelope,
            lambda ledger: apply_family_metadata(self.envelope, ledger, self.user),
        )

    def next_review(self):
        self.envelope.update(
            review_sequence=2, expected_previous_revision=self.envelope['revision']
        )
        self.envelope['individuals'][0]['metadata']['hpo_present'] = []
        seal(self.envelope)

    def test_create_explicit_parents_and_idempotence(self):
        first = self.apply()
        self.assertEqual(first.pk, self.apply().pk)
        self.assertEqual(Project.objects.count(), 1)
        self.assertEqual(Family.objects.count(), 1)
        self.assertEqual(Individual.objects.count(), 3)
        child = Individual.objects.get(individual_id='SYNTH_proband')
        self.assertEqual(child.mother.individual_id, 'SYNTH_mother')
        self.assertEqual(child.father.individual_id, 'SYNTH_father')
        self.assertEqual(child.sex, 'U')
        self.assertEqual(child.mother.sex, 'U')
        self.assertFalse(first.payload['searchable'])

    def test_review_metadata_clear_without_duplication(self):
        self.apply()
        self.next_review()
        self.apply()
        self.assertEqual(
            Individual.objects.get(individual_id='SYNTH_father').features, []
        )
        self.assertEqual(Individual.objects.count(), 3)
        self.assertEqual(PreparationLedger.objects.get().review_sequence, 2)

    def test_rollback_entire_family_on_third_update_failure(self):
        from seqr.views.utils.json_to_orm_utils import update_individual_from_json

        calls = []

        def fail_third(*args, **kwargs):
            calls.append(args[0].pk)
            if len(calls) == 3:
                raise RuntimeError('Synthetic third-person failure')
            return update_individual_from_json(*args, **kwargs)

        with patch(
            'seqr.views.utils.json_to_orm_utils.update_individual_from_json',
            side_effect=fail_third,
        ):
            with self.assertRaises(RuntimeError):
                self.apply()
        for model in (
            Individual,
            Family,
            Project,
            PreparationLedger,
            PreparationOutbox,
        ):
            self.assertFalse(model.objects.exists())

    def test_reject_unowned_existing_project(self):
        Project.objects.create(
            name=self.envelope['family_id'], genome_version='38', created_by=self.user
        )
        with self.assertRaises(PreparationConflict):
            self.apply()
        self.assertFalse(PreparationLedger.objects.exists())
        self.assertFalse(Individual.objects.exists())

    def test_reject_manual_metadata_and_parent_edits(self):
        self.apply()
        child = Individual.objects.get(individual_id='SYNTH_proband')
        child.mother = None
        child.save()
        self.next_review()
        with self.assertRaises(PreparationConflict):
            self.apply()
        self.assertEqual(PreparationLedger.objects.get().review_sequence, 1)
        self.assertEqual(PreparationOutbox.objects.count(), 1)

    def test_reject_loaded_active_or_inactive_data(self):
        self.apply()
        dataset = Dataset.objects.create(
            sample_type='WES',
            dataset_type=Dataset.DATASET_TYPE_VARIANT_CALLS,
            data_source='synthetic',
            loaded_date=timezone.now(),
            created_by=self.user,
        )
        person = Individual.objects.get(individual_id='SYNTH_mother')
        self.next_review()
        for relation in (dataset.active_individuals, dataset.inactive_individuals):
            relation.add(person)
            with self.assertRaises(PreparationConflict):
                self.apply()
            relation.clear()
        self.assertEqual(PreparationLedger.objects.get().review_sequence, 1)

    def test_reject_extra_family(self):
        self.apply()
        Family.objects.create(
            project=Project.objects.get(), family_id='UNOWNED', created_by=self.user
        )
        self.next_review()
        with self.assertRaises(PreparationConflict):
            self.apply()

    def test_reject_source_or_assay_replacement(self):
        self.apply()
        self.next_review()
        self.envelope['sample_type'] = 'WGS'
        self.envelope['request_id'] = canonical_digest(family_request(self.envelope))
        seal(self.envelope)
        with self.assertRaises(PreparationConflict):
            self.apply()

    def test_rollback_revision_preserves_previous_family(self):
        self.apply()
        previous = copy.deepcopy(
            list(Individual.objects.values('individual_id', 'features', 'sex'))
        )
        self.next_review()

        def update_then_fail(ledger):
            apply_family_metadata(self.envelope, ledger, self.user)
            raise RuntimeError('Synthetic post-update failure')

        with self.assertRaises(RuntimeError):
            apply_preparation(self.envelope, update_then_fail)
        self.assertEqual(
            list(Individual.objects.values('individual_id', 'features', 'sex')),
            previous,
        )
        self.assertEqual(PreparationLedger.objects.get().review_sequence, 1)

    @override_settings(PREPARATION_FAMILY_ENABLED=False)
    def test_disabled_by_default_gate(self):
        with self.assertRaises(PreparationConflict):
            self.apply()
        self.assertFalse(Project.objects.exists())
        self.assertFalse(PreparationLedger.objects.exists())
