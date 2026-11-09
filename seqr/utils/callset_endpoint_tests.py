import json
import socket
from types import SimpleNamespace
from unittest.mock import patch
from django.test import SimpleTestCase, RequestFactory, override_settings
from django.conf import settings
from seqr.views.apis import data_manager_api as api
from seqr.utils import callset_binding_tests as fixtures

ALIAS, PREPARED, ROOT = fixtures.ALIAS, fixtures.PREPARED, fixtures.ROOT


@override_settings(
    CALLSET_BINDING_ENABLED=True,
    CALLSET_LOADING_MAINTENANCE=False,
    LOADING_DATASETS_DIR=ROOT,
    CALLSET_PREPARED_PREFIX="lossless-v1",
    CALLSET_VALIDATION_MAX_AGE=3600,
)
class RealRequestTests(SimpleTestCase):
    def setUp(self):
        f = fixtures.CallsetBindingTests()
        f.setUp()
        self.fixture = f
        self.user = SimpleNamespace(
            pk=1,
            email="synthetic@example.invalid",
            is_authenticated=True,
            is_active=True,
            is_staff=True,
            is_superuser=True,
            userpolicy=SimpleNamespace(
                privacy_version=settings.SEQR_PRIVACY_VERSION,
                tos_version=settings.SEQR_TOS_VERSION,
            ),
        )
        self.body = {**f.request, "projects": []}
        self.patches = [
            patch.object(
                socket.socket,
                "connect",
                side_effect=AssertionError(
                    "Network forbidden in synthetic endpoint tests"
                ),
            ),
            patch("seqr.utils.callset_binding.storage.Client", return_value=f.client),
            patch.object(api, "LOADING_DATASETS_DIR", ROOT),
            patch.object(
                api.AirtableSession, "is_airtable_enabled", return_value=False
            ),
            patch.object(api, "anvil_enabled", return_value=False),
        ]
        for p in self.patches:
            p.start()
            self.addCleanup(p.stop)

    def request(self, body):
        r = RequestFactory().post(
            "/synthetic/", data=json.dumps(body), content_type="application/json"
        )
        r.user = self.user
        return r

    def validate(self):
        with patch.object(
            api, "validate_vcf_and_get_samples", return_value=["HEADER_SAMPLE"]
        ) as header:
            response = api.validate_callset(self.request(self.body))
        self.assertEqual(response.status_code, 200)
        self.assertEqual(header.call_args.args[0], PREPARED)
        return json.loads(response.content)

    def test_real_request_signing_alias_replacement_and_enqueue(self):
        result = self.validate()
        self.body.update(result)
        self.fixture.blobs[ALIAS].generation = "99"
        self.fixture.blobs[ALIAS].metadata = {}
        with (
            patch.object(api.Project.objects, "filter", return_value=[]),
            patch.object(api, "trigger_data_loading", return_value=True) as enqueue,
        ):
            response = api.load_data(self.request(self.body))
        self.assertEqual(response.status_code, 200)
        self.assertEqual(enqueue.call_args.args[5], PREPARED)
        self.assertNotIn("callset_identity", enqueue.call_args.kwargs)

    def test_actual_decorator_rejects_unauthorized_user(self):
        from django.core.exceptions import PermissionDenied

        self.user.is_staff = False
        self.user.is_superuser = False
        with self.assertRaises(PermissionDenied):
            api.validate_callset(self.request(self.body))

    def test_signed_payload_context_tamper_stops_before_database(self):
        self.body.update(self.validate())
        self.body["genomeVersion"] = "37"
        with (
            patch.object(api.Project.objects, "filter") as query,
            patch.object(api, "trigger_data_loading") as enqueue,
        ):
            self.assertEqual(api.load_data(self.request(self.body)).status_code, 400)
        query.assert_not_called()
        enqueue.assert_not_called()

    def test_generation_change_during_header_validation_never_signs(self):
        def header(*a, **k):
            self.fixture.blobs[PREPARED].generation = "999"
            return ["HEADER_SAMPLE"]

        with (
            patch.object(api, "validate_vcf_and_get_samples", side_effect=header),
            patch.object(api, "sign_validation") as sign,
        ):
            self.assertEqual(
                api.validate_callset(self.request(self.body)).status_code, 400
            )
        sign.assert_not_called()

    @override_settings(CALLSET_LOADING_MAINTENANCE=True)
    def test_maintenance_rejects_before_storage_database_or_enqueue(self):
        with (
            patch.object(api, "validated_callset") as binding,
            patch.object(api.Project.objects, "filter") as query,
            patch.object(api, "trigger_data_loading") as enqueue,
        ):
            self.assertEqual(api.load_data(self.request({})).status_code, 503)
        binding.assert_not_called()
        query.assert_not_called()
        enqueue.assert_not_called()

    def test_legacy_root_snapshot_survives_source_replacement(self):
        self.fixture.blobs[ALIAS].metadata = {}
        with patch.object(api, "validate_vcf_and_get_samples", return_value=["HEADER_SAMPLE"]) as header:
            response = api.validate_callset(self.request(self.body))
        self.assertEqual(response.status_code, 200)
        snapshot = header.call_args.args[0]
        self.assertTrue(snapshot.startswith(ROOT + ".seqr-snapshots/"))
        self.assertNotEqual(snapshot, ALIAS)
        self.assertEqual(len(self.fixture.rewrites), 1)
        self.body.update(json.loads(response.content))
        self.fixture.blobs[ALIAS].generation = "99"
        with (
            patch.object(api.Project.objects, "filter", return_value=[]),
            patch.object(api, "trigger_data_loading", return_value=True) as enqueue,
        ):
            response = api.load_data(self.request(self.body))
        self.assertEqual(response.status_code, 200)
        self.assertEqual(enqueue.call_args.args[5], snapshot)
        self.assertNotIn("callset_identity", enqueue.call_args.kwargs)
        self.assertEqual(len(self.fixture.rewrites), 1)

    def test_root_wildcard_never_reads_header_or_lists_objects(self):
        self.body["filePath"] = "/*.vcf.gz"
        with patch.object(api, "validate_vcf_and_get_samples") as header:
            response = api.validate_callset(self.request(self.body))
        self.assertEqual(response.status_code, 400)
        header.assert_not_called()
        self.fixture.client.list_blobs.assert_not_called()

    def test_dropdown_keeps_root_aliases_separate_and_hides_snapshots(self):
        from seqr.utils import vcf_utils

        with patch.object(vcf_utils, "list_files", return_value=[
            ROOT + "CASE_ONE.gatk.vcf.gz", ROOT + "CASE_TWO.gatk.vcf.gz",
            ROOT + ".seqr-snapshots/unique/input.vcf.gz",
        ]):
            self.assertEqual(vcf_utils.get_vcf_list(ROOT, self.user, preserve_root_aliases=True),
                             ["CASE_ONE.gatk.vcf.gz", "CASE_TWO.gatk.vcf.gz"])
