import json
from types import SimpleNamespace
from django.test import SimpleTestCase
from unittest.mock import Mock, patch


from django.test import override_settings

from seqr.utils.callset_binding import resolve_callset, sign_validation, validated_callset, verify_identity

ROOT = 'gs://synthetic-loading/'
ALIAS = ROOT + 'CASE_TEST.gatk.vcf.gz'
PREPARED = ROOT + 'lossless-v1/CASE_TEST/10-20/12345678-1234-1234-1234-123456789abc.gatk.vcf.gz'


@override_settings(CALLSET_BINDING_ENABLED=True, LOADING_DATASETS_DIR=ROOT,
                   CALLSET_PREPARED_PREFIX='lossless-v1', CALLSET_VALIDATION_MAX_AGE=3600)
class CallsetBindingTests(SimpleTestCase):
    def setUp(self):
        self.identity = {
            'vcf': {'path': PREPARED, 'generation': '30', 'size': '100', 'crc32c': 'AAAAAA==', 'sha256': 'a' * 64},
            'index': {'path': PREPARED + '.tbi', 'generation': '31', 'size': '20', 'crc32c': 'BAAAAA==', 'sha256': 'b' * 64},
        }
        binding = {'schema_version': 1, 'accession_id': 'CASE_TEST', **self.identity}
        self.blobs = {i['path']: SimpleNamespace(generation=i['generation'], size=int(i['size']), crc32c=i['crc32c'], metadata={}) for i in self.identity.values()}
        self.blobs[ALIAS] = SimpleNamespace(generation='40', size=100, crc32c='AAAAAA==', metadata={'seqr_binding': json.dumps(binding)})
        self.client = Mock()
        self.client.list_blobs.return_value = []
        self.client.bucket.return_value.get_blob.side_effect = lambda name: self.blobs.get(ROOT + name)
        self.copy_mode = 'success'
        self.rewrites = []
        self.client.bucket.return_value.blob.side_effect = self.copy_blob
        self.user = SimpleNamespace(pk=1)
        self.request = {'filePath': '/CASE_TEST.gatk.vcf.gz', 'genomeVersion': '38', 'sampleType': 'WES', 'vcfSamples': ['HEADER_SAMPLE']}

    def copy_blob(self, name, generation=None):
        if generation is not None:
            return SimpleNamespace(name=name, generation=generation)
        target = SimpleNamespace(metadata=None)
        def rewrite(source, **kwargs):
            self.rewrites.append((name, source, kwargs))
            self.assertEqual(kwargs['if_generation_match'], 0)
            self.assertEqual(kwargs['if_source_generation_match'], source.generation)
            self.assertIsNone(kwargs['retry'])
            if self.copy_mode == 'interrupt' or (self.copy_mode == 'index_interrupt' and name.endswith('.tbi')):
                raise TimeoutError('synthetic interruption')
            original = self.blobs[ROOT + source.name]
            if self.copy_mode == 'source_changed':
                from google.api_core.exceptions import PreconditionFailed
                original.generation = '999'
                raise PreconditionFailed('Pinned source generation no longer available')
            self.assertEqual(str(source.generation), str(original.generation))
            if self.copy_mode == 'continuation':
                if kwargs['token'] is None:
                    return 'token-one', 20, original.size
                if kwargs['token'] == 'token-one':
                    return 'token-two', 70, original.size
                self.assertEqual(kwargs['token'], 'token-two')
            self.blobs[ROOT + name] = SimpleNamespace(generation='101', size=original.size,
                crc32c=original.crc32c, metadata=target.metadata)
            if self.copy_mode.startswith('mismatch_'):
                completed = self.blobs[ROOT + name]
                field = self.copy_mode.removeprefix('mismatch_')
                setattr(completed, field, {'metadata': {}, 'size': 999, 'crc32c': 'BAAAAA=='}[field])
                raise TimeoutError('synthetic lost response with mismatched destination')
            if self.copy_mode == 'lost_response':
                raise TimeoutError('synthetic lost response')
            return None, original.size, original.size
        target.rewrite = rewrite
        return target

    def token(self):
        self.request['callsetValidation'] = sign_validation(self.request, self.user, ALIAS, self.identity, ['HEADER_SAMPLE'])

    def test_alias_resolution_and_repeated_request(self):
        self.assertEqual(resolve_callset(ALIAS, self.client), (PREPARED, self.identity))
        self.token()
        for _ in range(2):
            self.assertEqual(validated_callset(self.request, self.user, ALIAS, self.client), (PREPARED, self.identity))

    def test_alias_replaced_after_validation_does_not_change_selected_bytes(self):
        self.token()
        self.blobs[ALIAS].generation = '99'
        self.blobs[ALIAS].metadata = {}
        self.assertEqual(validated_callset(self.request, self.user, ALIAS, self.client)[0], PREPARED)

    def test_same_size_replacement_fails_for_either_object(self):
        self.token()
        for name in ('vcf', 'index'):
            blob = self.blobs[self.identity[name]['path']]
            original = blob.generation
            blob.generation = '99'
            with self.assertRaises(ValueError):
                validated_callset(self.request, self.user, ALIAS, self.client)
            blob.generation = original

    def test_forged_missing_expired_token(self):
        for token in ('', 'forged'):
            self.request['callsetValidation'] = token
            with self.assertRaises(ValueError):
                validated_callset(self.request, self.user, ALIAS, self.client)
        with patch('django.core.signing.time.time', return_value=1):
            self.token()
        with self.assertRaises(ValueError):
            validated_callset(self.request, self.user, ALIAS, self.client)

    def test_context_mismatch(self):
        self.token()
        for key, value in [('sampleType', 'WGS'), ('genomeVersion', '37'), ('datasetType', 'SV'), ('vcfSamples', ['OTHER'])]:
            request = {**self.request, key: value}
            with self.assertRaises(ValueError):
                validated_callset(request, self.user, ALIAS, self.client)
        for user, path in [(SimpleNamespace(pk=2), ALIAS), (self.user, ROOT + 'OTHER.gatk.vcf.gz')]:
            with self.assertRaises(ValueError):
                validated_callset(self.request, user, path, self.client)

    def test_missing_malformed_and_cross_bucket_alias_bindings(self):
        for metadata in ({'seqr_binding': None}, {'seqr_binding': '{'}, {'seqr_binding': 'null'}, {'seqr_binding': '{}'}):
            self.blobs[ALIAS].metadata = metadata
            with self.assertRaises(ValueError):
                resolve_callset(ALIAS, self.client)
        identity = {**self.identity, 'vcf': {**self.identity['vcf'], 'path': 'gs://untrusted/file.vcf.gz'}}
        self.blobs[ALIAS].metadata = {'seqr_binding': json.dumps({'schema_version': 1, 'accession_id': 'CASE_TEST', **identity})}
        with self.assertRaises(ValueError):
            resolve_callset(ALIAS, self.client)

    def test_missing_pair_and_alias_byte_mismatch(self):
        self.blobs[ALIAS].crc32c = 'changed'
        with self.assertRaises(ValueError):
            resolve_callset(ALIAS, self.client)
        del self.blobs[PREPARED + '.tbi']
        with self.assertRaises(ValueError):
            verify_identity(self.identity, self.client)

    def test_unique_legacy_path_is_bound(self):
        path, identity = resolve_callset(PREPARED, self.client)
        self.assertEqual(path, PREPARED)
        self.assertEqual(identity['vcf']['generation'], '30')
        self.assertNotIn('sha256', identity['vcf'])

    def test_old_root_gatk_without_metadata_is_generation_bound(self):
        self.blobs[ALIAS].metadata = {}
        resolved, identity = resolve_callset(ALIAS, self.client)
        self.assertIn('/.seqr-snapshots/', resolved)
        self.assertEqual(identity['vcf']['generation'], '101')
        self.assertIsNone(identity['index'])

    def test_legacy_root_vcf_is_bound_without_index(self):
        path = ROOT + 'legacy.vcf'
        self.blobs[path] = SimpleNamespace(generation='71', size=50, crc32c='AAAAAA==', metadata={})
        resolved, identity = resolve_callset(path, self.client)
        self.assertIn('/.seqr-snapshots/', resolved)
        self.assertIsNone(identity['index'])
        self.request['callsetValidation'] = sign_validation(self.request, self.user, path, identity, ['HEADER_SAMPLE'])
        self.assertEqual(validated_callset(self.request, self.user, path, self.client), (resolved, identity))
        self.blobs[path].generation = '72'
        self.assertEqual(validated_callset(self.request, self.user, path, self.client), (resolved, identity))

    def test_source_generation_changes_between_precheck_and_rewrite(self):
        self.blobs[ALIAS].metadata = {}
        self.copy_mode = 'source_changed'
        with self.assertRaises(ValueError):
            resolve_callset(ALIAS, self.client)
        self.assertEqual(len(self.rewrites), 1)
        self.assertEqual(self.rewrites[0][2]['if_source_generation_match'], 40)
        self.assertEqual(self.blobs[ALIAS].generation, '999')
        self.assertFalse(any('/.seqr-snapshots/' in path for path in self.blobs))

    def test_multitoken_rewrite_retains_pins_and_never_retries(self):
        self.blobs[ALIAS].metadata = {}
        self.copy_mode = 'continuation'
        path, identity = resolve_callset(ALIAS, self.client)
        self.assertIn('/.seqr-snapshots/', path)
        self.assertEqual(identity['vcf']['generation'], '101')
        self.assertEqual([kwargs['token'] for _, _, kwargs in self.rewrites], [None, 'token-one', 'token-two'])
        for name, source, kwargs in self.rewrites:
            self.assertEqual(ROOT + name, path)
            self.assertEqual(source.generation, 40)
            self.assertEqual(kwargs['if_source_generation_match'], 40)
            self.assertEqual(kwargs['if_generation_match'], 0)
            self.assertIsNone(kwargs['retry'])

    def test_lost_response_reconciliation_rejects_mismatched_marker_or_bytes(self):
        self.blobs[ALIAS].metadata = {}
        for field in ('metadata', 'size', 'crc32c'):
            with self.subTest(field=field):
                self.copy_mode = 'mismatch_' + field
                self.rewrites.clear()
                with self.assertRaises(ValueError):
                    resolve_callset(ALIAS, self.client)
                self.assertEqual(len(self.rewrites), 1)

    def test_snapshot_interruption_never_returns_unbound_path(self):
        self.blobs[ALIAS].metadata = {}
        self.copy_mode = 'interrupt'
        with self.assertRaises(ValueError):
            resolve_callset(ALIAS, self.client)
        self.assertEqual(len(self.rewrites), 1)

    def test_partial_pair_snapshot_never_returns_or_retries(self):
        self.blobs[ALIAS].metadata = {}
        self.blobs[ALIAS + '.tbi'] = SimpleNamespace(generation='41', size=20, crc32c='BAAAAA==', metadata={})
        self.copy_mode = 'index_interrupt'
        with self.assertRaises(ValueError):
            resolve_callset(ALIAS, self.client)
        self.assertEqual(len(self.rewrites), 2)
        self.assertEqual(len([path for path in self.blobs if '/.seqr-snapshots/' in path]), 1)

    def test_snapshot_lost_response_reconciles_without_repeating_write(self):
        self.blobs[ALIAS].metadata = {}
        self.copy_mode = 'lost_response'
        path, identity = resolve_callset(ALIAS, self.client)
        self.assertIn('/.seqr-snapshots/', path)
        self.assertEqual(identity['vcf']['generation'], '101')
        self.assertEqual(len(self.rewrites), 1)

    def test_snapshot_pair_and_revalidation_use_fresh_paths(self):
        self.blobs[ALIAS].metadata = {}
        self.blobs[ALIAS + '.tbi'] = SimpleNamespace(generation='41', size=20, crc32c='BAAAAA==', metadata={})
        path, identity = resolve_callset(ALIAS, self.client)
        self.assertEqual(identity['index']['path'], path + '.tbi')
        self.blobs[ALIAS].generation = '42'  # equal size, new content/version
        next_path, _ = resolve_callset(ALIAS, self.client)
        self.assertNotEqual(path, next_path)
        self.assertEqual(len(self.rewrites), 4)

    def test_empty_root_glob_rejects_future_aliases(self):
        for path in (ROOT + '*.vcf.gz', ROOT + 'CASE*.gatk.vcf.gz', ROOT + '**/*.vcf.gz'):
            with self.assertRaises(ValueError):
                resolve_callset(path, self.client)
        self.client.list_blobs.assert_not_called()

    def test_unmanaged_wildcard_and_matrix_table_remain_legacy(self):
        for path in (ROOT + 'shards/*.vcf.gz', ROOT + 'callset.mt', 'gs://other/legacy.vcf.gz'):
            self.assertEqual(resolve_callset(path, self.client), (path, None))
            self.assertEqual(validated_callset({}, self.user, path, self.client), (path, None))
        self.client.bucket.assert_not_called()

    def test_broad_legacy_glob_cannot_include_managed_alias(self):
        self.client.list_blobs.return_value = [SimpleNamespace(name='CASE_TEST.gatk.vcf.gz', metadata={'seqr_binding': None})]
        with self.assertRaises(ValueError):
            resolve_callset(ROOT + '*.vcf.gz', self.client)

    def test_unsafe_paths_rejected(self):
        for path in (ROOT + '../outside.vcf.gz', ROOT + 'name;cmd.vcf.gz', ROOT + '/double.vcf.gz'):
            with self.assertRaises(ValueError):
                resolve_callset(path, self.client)

    @override_settings(CALLSET_BINDING_ENABLED=False)
    def test_disabled_legacy_flow_does_not_access_storage(self):
        self.assertEqual(resolve_callset(ALIAS, self.client), (ALIAS, None))
        self.assertEqual(validated_callset({}, self.user, ALIAS, self.client), (ALIAS, None))
        self.client.bucket.assert_not_called()


class CallsetEndpointContractTests(SimpleTestCase):
    """Execute endpoint bodies without a database or production credentials."""
    def endpoint(self, name, **namespace):
        import ast
        from pathlib import Path
        source = Path(__file__).parents[1] / 'views/apis/data_manager_api.py'
        tree = ast.parse(source.read_text())
        function = next(n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == name)
        function.decorator_list = []
        context = {'settings': SimpleNamespace(), 'json': json, 'create_json_response': lambda value, status=200: (status, value), **namespace}
        exec(compile(ast.Module(body=[function], type_ignores=[]), str(source), 'exec'), context)
        return context[name]

    def test_validation_returns_token_and_checks_version_after_header(self):
        order = []
        verify = Mock(side_effect=lambda identity: order.append('verify'))
        validate = Mock(side_effect=lambda *args, **kwargs: order.append('header') or {'HEADER_SAMPLE'})
        endpoint = self.endpoint('validate_callset', anvil_enabled=lambda: False,
            _callset_path=lambda request: ALIAS, enabled_for=lambda path: True,
            resolve_callset=lambda path: (PREPARED, {'identity': 'test'}),
            validate_vcf_and_get_samples=validate, verify_identity=verify,
            sign_validation=lambda *args: 'signed-token')
        status, response = endpoint(SimpleNamespace(body=json.dumps({'genomeVersion': '38', 'filePath': ALIAS}), user=object()))
        self.assertEqual(status, 200)
        self.assertEqual(response['callsetValidation'], 'signed-token')
        self.assertEqual(order, ['header', 'verify'])
        self.assertEqual(validate.call_args.args[0], PREPARED)

    def test_changed_generation_after_header_never_signs(self):
        sign = Mock()
        endpoint = self.endpoint('validate_callset', anvil_enabled=lambda: False,
            _callset_path=lambda request: ALIAS, enabled_for=lambda path: True,
            resolve_callset=lambda path: (PREPARED, {'identity': 'test'}),
            validate_vcf_and_get_samples=lambda *args, **kwargs: {'HEADER_SAMPLE'},
            verify_identity=Mock(side_effect=ValueError('changed')), sign_validation=sign)
        status, _ = endpoint(SimpleNamespace(body=json.dumps({'genomeVersion': '38', 'filePath': ALIAS}), user=object()))
        self.assertEqual(status, 400)
        sign.assert_not_called()

    def test_load_rejects_invalid_token_before_project_lookup_or_enqueue(self):
        trigger = Mock()
        endpoint = self.endpoint('load_data', _callset_path=lambda request: ALIAS,
            validated_callset=Mock(side_effect=ValueError('expired')), trigger_data_loading=trigger)
        status, _ = endpoint(SimpleNamespace(body='{}', user=object()))
        self.assertEqual(status, 400)
        trigger.assert_not_called()
