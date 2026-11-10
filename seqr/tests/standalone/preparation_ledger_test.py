"""Run standalone only: SQLite transaction tests, not PostgreSQL lock certification."""
import copy
import hashlib
import json
import os
import unittest

from django.conf import settings
database = {'ENGINE': 'django.db.backends.sqlite3', 'NAME': ':memory:'}
if os.environ.get('PREP_TEST_POSTGRES_PORT'):
    database = {'ENGINE': 'django.db.backends.postgresql', 'NAME': 'prep_ledger_synthetic', 'USER': 'synthetic', 'PASSWORD': 'synthetic', 'HOST': '127.0.0.1', 'PORT': int(os.environ['PREP_TEST_POSTGRES_PORT'])}
settings.configure(INSTALLED_APPS=[], DATABASES={'default': database}, USE_TZ=True)
import django
django.setup()
from django.db import connection, models
from seqr.preparation_models import PreparationLedger, PreparationOutbox
from seqr.utils.preparation_ledger import apply_preparation, canonical_digest, deliver_acknowledgment, PreparationConflict


class SyntheticClinicalRow(models.Model):
    value = models.CharField(max_length=64)
    class Meta:
        app_label = 'seqr'


class LedgerTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        with connection.schema_editor() as editor:
            for model in (PreparationLedger, PreparationOutbox, SyntheticClinicalRow):
                editor.create_model(model)

    def setUp(self):
        PreparationOutbox.objects.all().delete()
        PreparationLedger.objects.all().delete()
        SyntheticClinicalRow.objects.all().delete()
        self.envelope = dict(schema_version=2, kind='individual', accession_id='SYNTH1',
                             request_id='a' * 64, review_sequence=1, expected_previous_revision=None,
                             prepared=False, searchable=False, manual_load_required=True, genome_version='38', sample_type='WES',
                             individual=dict(sample_id='SYNTH1', sex='U', affected='U', metadata=dict(hpo_present=[], hpo_absent=[])))
        self.envelope.update(source=dict(bucket='synthetic', path='fixture.vcf', generation='1', vcf_sample_id='OLD'), source_sha256='c' * 64, records=1, status='awaiting_consumer')
        request = {key: self.envelope[key] for key in ('accession_id', 'sample_type', 'genome_version')}
        request = {'schema_version': 1, 'kind': 'normalize_singleton', **request, 'source': self.envelope['source']}
        self.envelope['request_id'] = hashlib.sha256(json.dumps(request, separators=(',', ':')).encode()).hexdigest()
        self.seal()
        self.calls = 0

    def seal(self):
        self.envelope['revision'] = canonical_digest({k: v for k, v in self.envelope.items() if k != 'revision'})

    def apply(self, ledger):
        self.calls += 1
        SyntheticClinicalRow.objects.create(value=self.envelope['revision'])
        return 'R_SYNTH', 'b' * 64

    def test_singleton_cannot_claim_reserved_family_project_name(self):
        from seqr.utils.preparation_family_test_utils import synthetic_family_envelope
        reserved = synthetic_family_envelope()['family_id']
        self.envelope['accession_id'] = reserved
        self.envelope['individual']['sample_id'] = reserved
        request = {'schema_version': 1, 'kind': 'normalize_singleton',
                   **{key: self.envelope[key] for key in ('accession_id', 'sample_type', 'genome_version')},
                   'source': self.envelope['source']}
        self.envelope['request_id'] = hashlib.sha256(json.dumps(request, separators=(',', ':')).encode()).hexdigest()
        self.seal()
        with self.assertRaisesRegex(PreparationConflict, 'Reserved family'):
            apply_preparation(self.envelope, self.apply)
        self.assertFalse(PreparationLedger.objects.exists())

    def test_family_namespace_ack_and_metadata_revision(self):
        from seqr.utils.preparation_family_test_utils import synthetic_family_envelope, seal
        family = synthetic_family_envelope()
        singleton = apply_preparation(self.envelope, self.apply)
        first = apply_preparation(family, self.apply)
        self.assertNotEqual(first.ledger_id, singleton.ledger_id)
        self.assertEqual(first.ledger.subject_key, family['subject_key'])
        self.assertEqual(first.payload['joint_report_id'], family['joint_report_id'])
        self.assertEqual(first.payload['subject_key'], family['subject_key'])
        self.assertEqual(first.payload['kind'], 'family_joint_callset')
        self.assertNotIn('accession_id', first.payload)
        self.assertEqual(apply_preparation(family, self.apply).pk, first.pk)
        family.update(review_sequence=2, expected_previous_revision=family['revision'])
        family['individuals'][0]['metadata']['hpo_present'] = []
        seal(family)
        second = apply_preparation(family, self.apply)
        self.assertEqual(second.review_sequence, 2)
        self.assertEqual(self.calls, 3)

    def test_family_source_change_and_callback_rollback(self):
        from seqr.utils.preparation_family_test_utils import synthetic_family_envelope, seal
        from seqr.utils.preparation_family_contract import family_request
        family = synthetic_family_envelope()
        def fail(ledger):
            self.apply(ledger)
            raise RuntimeError('synthetic family callback failure')
        with self.assertRaises(RuntimeError):
            apply_preparation(family, fail)
        self.assertFalse(PreparationLedger.objects.exists())
        self.assertFalse(SyntheticClinicalRow.objects.exists())
        apply_preparation(family, self.apply)
        family.update(review_sequence=2, expected_previous_revision=family['revision'])
        family['source']['generation'] = '2'
        family['request_id'] = canonical_digest(family_request(family))
        seal(family)
        with self.assertRaises(PreparationConflict):
            apply_preparation(family, self.apply)
        self.assertEqual(PreparationLedger.objects.get().review_sequence, 1)

    def test_apply_and_exact_retry_once(self):
        first = apply_preparation(self.envelope, self.apply)
        second = apply_preparation(self.envelope, self.apply)
        self.assertEqual(first.pk, second.pk)
        self.assertEqual(self.calls, 1)
        self.assertEqual(SyntheticClinicalRow.objects.count(), 1)
        self.assertFalse(first.payload['searchable'])

    def test_failure_rolls_back_clinical_and_ledger(self):
        def fail(ledger):
            self.apply(ledger)
            raise RuntimeError('synthetic failure')
        with self.assertRaises(RuntimeError): apply_preparation(self.envelope, fail)
        self.assertEqual(SyntheticClinicalRow.objects.count(), 0)
        self.assertEqual(PreparationLedger.objects.count(), 0)
        self.assertEqual(PreparationOutbox.objects.count(), 0)

    def test_ack_failure_replays_without_clinical_reapply(self):
        event = apply_preparation(self.envelope, self.apply)
        sent = []
        def publish_then_crash(payload):
            sent.append(payload)
            raise RuntimeError('crash after external publish')
        with self.assertRaises(RuntimeError): deliver_acknowledgment(event.pk, publish_then_crash)
        event.refresh_from_db()
        self.assertIsNone(event.delivered_at)
        deliver_acknowledgment(event.pk, sent.append)
        self.assertEqual(sent[0], sent[1])
        self.assertEqual(self.calls, 1)
        deliver_acknowledgment(event.pk, sent.append)
        self.assertEqual(len(sent), 2)

    def test_ordered_revision_and_stale_rejection(self):
        first = copy.deepcopy(self.envelope)
        apply_preparation(first, self.apply)
        self.envelope.update(review_sequence=2, expected_previous_revision=first['revision'])
        self.seal()
        apply_preparation(self.envelope, self.apply)
        old_ack = apply_preparation(first, self.apply)
        self.assertEqual(old_ack.review_sequence, 1)
        self.assertEqual(PreparationLedger.objects.get().review_sequence, 2)
        self.assertEqual(self.calls, 2)

    def test_conflicting_revision_or_gap_rejected(self):
        apply_preparation(self.envelope, self.apply)
        self.envelope['review_sequence'] = 3
        self.envelope['expected_previous_revision'] = self.envelope['revision']
        self.seal()
        with self.assertRaises(PreparationConflict): apply_preparation(self.envelope, self.apply)
        self.assertEqual(self.calls, 1)

    @unittest.skipUnless(connection.vendor == 'postgresql', 'Requires isolated PostgreSQL row-lock semantics')
    def test_concurrent_duplicate_is_applied_once(self):
        from concurrent.futures import ThreadPoolExecutor
        from threading import Barrier
        from django.db import connections
        barrier = Barrier(2)
        def attempt():
            try:
                barrier.wait(timeout=5)
                return apply_preparation(self.envelope, self.apply).pk
            finally:
                connections.close_all()
        with ThreadPoolExecutor(max_workers=2) as pool:
            first, second = list(pool.map(lambda _: attempt(), range(2)))
        self.assertEqual(first, second)
        self.assertEqual(self.calls, 1)
        self.assertEqual(SyntheticClinicalRow.objects.count(), 1)

    @unittest.skipUnless(connection.vendor == 'postgresql', 'Requires isolated PostgreSQL row-lock semantics')
    def test_concurrent_conflicting_successors_only_one_commits(self):
        from concurrent.futures import ThreadPoolExecutor
        from threading import Barrier
        from django.db import connections
        apply_preparation(self.envelope, self.apply)
        predecessor = self.envelope['revision']
        barrier = Barrier(2)
        candidates = []
        for sex in ('M', 'F'):
            candidate = copy.deepcopy(self.envelope)
            candidate.update(review_sequence=2, expected_previous_revision=predecessor)
            candidate['individual']['sex'] = sex
            candidate['revision'] = canonical_digest({k: v for k, v in candidate.items() if k != 'revision'})
            candidates.append(candidate)
        def attempt(candidate):
            try:
                barrier.wait(timeout=5)
                return apply_preparation(candidate, self.apply).revision
            except PreparationConflict:
                return None
            finally:
                connections.close_all()
        with ThreadPoolExecutor(max_workers=2) as pool:
            outcomes = list(pool.map(attempt, candidates))
        self.assertEqual(sum(value is not None for value in outcomes), 1)
        self.assertEqual(PreparationLedger.objects.get().review_sequence, 2)
        self.assertEqual(PreparationOutbox.objects.count(), 2)
        self.assertEqual(SyntheticClinicalRow.objects.count(), 2)

    def test_source_and_assay_remain_pinned_across_metadata_reviews(self):
        apply_preparation(self.envelope, self.apply)
        original = copy.deepcopy(self.envelope)
        for field, value in [('source', dict(original['source'], generation='2')), ('sample_type', 'WGS')]:
            self.envelope = copy.deepcopy(original)
            self.envelope.update(review_sequence=2, expected_previous_revision=original['revision'])
            self.envelope[field] = value
            request = {key: self.envelope[key] for key in ('accession_id', 'sample_type', 'genome_version')}
            request = {'schema_version': 1, 'kind': 'normalize_singleton', **request, 'source': self.envelope['source']}
            self.envelope['request_id'] = hashlib.sha256(json.dumps(request, separators=(',', ':')).encode()).hexdigest()
            self.seal()
            with self.subTest(field=field), self.assertRaises(PreparationConflict):
                apply_preparation(self.envelope, self.apply)
        self.assertEqual(self.calls, 1)
        self.assertEqual(PreparationLedger.objects.get().request_id, original['request_id'])

    def test_tampered_envelope_rejected(self):
        self.envelope['accession_id'] = 'OTHER'
        with self.assertRaises(PreparationConflict): apply_preparation(self.envelope, self.apply)
        self.assertEqual(self.calls, 0)


if __name__ == '__main__':
    unittest.main()
