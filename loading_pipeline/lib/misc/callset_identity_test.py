import copy
from unittest.mock import Mock, patch

import luigi
import luigi.util
import pytest

from loading_pipeline.lib.misc.callset_identity import (
    ensure_run_binding,
    guard_callset,
    identity_cache_key,
    normalize_identity,
    reject_unbound_managed_alias,
    verify_callset_identity,
)

PATH = 'gs://synthetic-inputs/prepared/synthetic.vcf.gz'
IDENTITY = {
    'vcf': {'path': PATH, 'generation': '101', 'size': '100', 'crc32c': 'AAAAAA=='},
    'index': {
        'path': PATH + '.tbi',
        'generation': '102',
        'size': '20',
        'crc32c': 'AAAAAA==',
    },
}


@pytest.fixture(autouse=True)
def disabled_binding(monkeypatch):
    monkeypatch.delenv('CALLSET_BINDING_ENABLED', raising=False)


def client_for(identity):
    client = Mock()

    def get(**kwargs):
        expected = next(
            value
            for value in identity.values()
            if value['path'] == 'gs://' + kwargs['bucket'] + '/' + kwargs['object']
        )
        assert kwargs['ifGenerationMatch'] == expected['generation']
        request = Mock()
        request.execute.return_value = expected
        return request

    client.objects.return_value.get.side_effect = get
    return client


def test_verify_exact_generations():
    client = client_for(IDENTITY)
    verify_callset_identity(PATH, IDENTITY, client)
    assert client.objects.return_value.get.call_count == len(IDENTITY)


@pytest.mark.parametrize('role', ['vcf', 'index'])
@pytest.mark.parametrize(
    'field,value',
    [('generation', '103'), ('crc32c', 'AAAAAQ=='), ('size', '101')],
)
def test_replacement_rejected(role, field, value):
    client = client_for(IDENTITY)
    changed = copy.deepcopy(IDENTITY)
    changed[role][field] = value
    responses = [changed['vcf'], changed['index']]
    client.objects.return_value.get.side_effect = None
    client.objects.return_value.get.return_value.execute.side_effect = responses
    with pytest.raises(ValueError, match='changed'):
        verify_callset_identity(PATH, IDENTITY, client)


def test_missing_index_fails_closed():
    client = client_for(IDENTITY)
    client.objects.return_value.get.side_effect = [
        Mock(execute=Mock(return_value=IDENTITY['vcf'])),
        FileNotFoundError(),
    ]
    with pytest.raises(FileNotFoundError):
        verify_callset_identity(PATH, IDENTITY, client)


def test_same_size_generation_changes_cache_key():
    changed = copy.deepcopy(IDENTITY)
    changed['vcf']['generation'] = '103'
    assert identity_cache_key(PATH, changed) != identity_cache_key(PATH, IDENTITY)
    changed = copy.deepcopy(IDENTITY)
    changed['index']['generation'] = '104'
    assert identity_cache_key(PATH, changed) != identity_cache_key(PATH, IDENTITY)


@pytest.mark.parametrize(
    'field,value',
    [
        ('generation', 101),
        ('generation', '0'),
        ('size', '-1'),
        ('size', 100),
        ('crc32c', 'invalid'),
        ('sha256', 'xyz'),
        ('path', 'gs://another/input.vcf.gz'),
    ],
)
def test_malformed_identity_rejected(field, value):
    invalid = copy.deepcopy(IDENTITY)
    invalid['vcf'][field] = value
    with pytest.raises(ValueError):
        normalize_identity(PATH, invalid)


@pytest.mark.parametrize(
    'identity',
    [{'vcf': {}}, {'vcf': None, 'index': {}}, ['invalid'], [], None],
)
def test_malformed_pair_rejected(identity):
    with pytest.raises(ValueError):
        normalize_identity(PATH, identity)


def test_legacy_compatible_and_flag_fail_closed(monkeypatch):
    verify_callset_identity('legacy.mt', {}, client=Mock())
    monkeypatch.setenv('CALLSET_BINDING_ENABLED', 'true')
    verify_callset_identity('legacy.mt', {})
    verify_callset_identity(PATH, IDENTITY, client_for(IDENTITY))


class Params(luigi.Task):
    callset_path = luigi.Parameter(default=PATH)
    callset_identity = luigi.DictParameter(default=IDENTITY)


@guard_callset
@luigi.util.inherits(Params)
class MaterializedTask(luigi.Task):
    def complete(self):
        return True

    def run(self):
        self.materialize()

    def materialize(self):
        pass


def test_luigi_frozen_identity_and_serialization():
    task = MaterializedTask()
    assert normalize_identity(task.callset_path, task.callset_identity) == IDENTITY
    reloaded = MaterializedTask.from_str_params(task.to_str_params())
    assert dict(reloaded.callset_identity) == IDENTITY


def test_completed_cache_checked_before_reuse():
    with (
        patch(
            'loading_pipeline.lib.misc.callset_identity.verify_callset_identity',
            side_effect=ValueError('replaced'),
        ),
        pytest.raises(ValueError, match='replaced'),
    ):
        MaterializedTask().complete()


def test_materialized_import_fenced_before_and_after():
    events = []
    with (
        patch(
            'loading_pipeline.lib.misc.callset_identity.verify_callset_identity',
            side_effect=lambda *_: events.append('verify'),
        ),
        patch.object(
            MaterializedTask,
            'materialize',
            side_effect=lambda: events.append('materialize'),
        ),
    ):
        MaterializedTask().run()
    assert events == ['verify', 'materialize', 'verify']


def test_change_during_import_prevents_success_and_resume():
    with (
        patch(
            'loading_pipeline.lib.misc.callset_identity.verify_callset_identity',
            side_effect=[None, ValueError('replaced')],
        ),
        pytest.raises(ValueError, match='replaced'),
    ):
        MaterializedTask().run()
    with (
        patch(
            'loading_pipeline.lib.misc.callset_identity.verify_callset_identity',
            side_effect=ValueError('replaced'),
        ),
        pytest.raises(ValueError, match='replaced'),
    ):
        MaterializedTask().complete()


def test_generator_resume_checked():
    @guard_callset
    class GeneratorTask:
        callset_path = PATH
        callset_identity = IDENTITY

        def complete(self):
            return False

        def run(self):
            yield 'dependency'
            pytest.fail('Must not resume after replacement')

    with patch(
        'loading_pipeline.lib.misc.callset_identity.verify_callset_identity',
        side_effect=[None, None, ValueError('replaced')],
    ):
        gen = GeneratorTask().run()
        assert next(gen) == 'dependency'
        with pytest.raises(ValueError, match='replaced'):
            next(gen)


def test_generator_passes_luigi_dependency_result():
    received = []

    @guard_callset
    class GeneratorTask:
        callset_path = PATH
        callset_identity = IDENTITY

        def complete(self):
            return False

        def run(self):
            received.append((yield 'dependency'))

    with patch('loading_pipeline.lib.misc.callset_identity.verify_callset_identity'):
        generator = GeneratorTask().run()
        assert next(generator) == 'dependency'
        with pytest.raises(StopIteration):
            generator.send('completed target')
    assert received == ['completed target']


def test_legacy_vcf_without_index():
    identity = {'vcf': IDENTITY['vcf'], 'index': None}
    client = client_for(IDENTITY)
    verify_callset_identity(PATH, identity, client)
    assert client.objects.return_value.get.call_count == 1
    assert identity_cache_key(PATH, identity) != identity_cache_key(PATH, IDENTITY)


@pytest.mark.parametrize('binding', ['', '{}', None, 'malformed'])
def test_unbound_managed_metadata_rejected(monkeypatch, binding):
    monkeypatch.setenv('CALLSET_BINDING_ENABLED', 'true')
    client = Mock()
    client.objects.return_value.get.return_value.execute.return_value = {
        'metadata': {'seqr_binding': binding},
    }
    with pytest.raises(ValueError, match='Managed'):
        verify_callset_identity(PATH, {}, client)


def test_unbound_legacy_vcf_still_works(monkeypatch):
    monkeypatch.setenv('CALLSET_BINDING_ENABLED', 'true')
    client = Mock()
    client.objects.return_value.get.return_value.execute.return_value = {'metadata': {}}
    verify_callset_identity(PATH, {}, client)


def test_managed_wildcard_member_rejected(monkeypatch):
    import sys
    from types import ModuleType, SimpleNamespace

    monkeypatch.setenv('CALLSET_BINDING_ENABLED', 'true')
    fs = ModuleType('hailtop.fs')
    fs.ls = Mock(return_value=[SimpleNamespace(path=PATH)])
    hailtop = ModuleType('hailtop')
    hailtop.fs = fs
    client = Mock()
    client.objects.return_value.get.return_value.execute.return_value = {
        'metadata': {'seqr_binding': '{}'},
    }
    with patch.dict(sys.modules, {'hailtop': hailtop, 'hailtop.fs': fs}):
        with pytest.raises(ValueError, match='Managed'):
            reject_unbound_managed_alias('gs://synthetic-inputs/*.vcf.gz', client)
        client.objects.return_value.get.return_value.execute.return_value = {}
        reject_unbound_managed_alias('gs://synthetic-inputs/*.vcf.gz', client)


def test_local_run_binding_retries_and_collision(tmp_path):
    run = str(tmp_path / 'synthetic-run')
    ensure_run_binding(run, {'identity': IDENTITY, 'sample_type': 'WGS'})
    ensure_run_binding(run, {'identity': IDENTITY, 'sample_type': 'WGS'})
    with pytest.raises(ValueError, match='different'):
        ensure_run_binding(run, {'identity': IDENTITY, 'sample_type': 'WES'})
    changed = copy.deepcopy(IDENTITY)
    changed['vcf']['generation'] = '103'
    with pytest.raises(ValueError, match='different'):
        ensure_run_binding(run, {'identity': changed, 'sample_type': 'WGS'})


def test_existing_unbound_run_not_blessed(tmp_path):
    (tmp_path / '_SUCCESS').write_text('')
    with pytest.raises(ValueError, match='unbound'):
        ensure_run_binding(str(tmp_path), {'identity': IDENTITY})
    assert not (tmp_path / '_CALLSET_BINDING.json').exists()


class FakeGcs:
    def __init__(self, old_outputs=False, lost_write=False, collision=None):
        self.manifest = None
        self.old_outputs = old_outputs
        self.lost_write = lost_write
        self.collision = collision
        self.insert_calls = 0

    def objects(self):
        return self

    def get(self, **kwargs):
        def execute():
            if self.manifest is None:
                from googleapiclient.errors import HttpError
                from httplib2 import Response

                raise HttpError(Response({'status': '404'}), b'not found')
            return self.manifest

        return Mock(execute=execute)

    def list(self, **kwargs):
        return Mock(
            execute=Mock(return_value={'items': [{}]} if self.old_outputs else {}),
        )

    def insert(self, **kwargs):
        assert kwargs['ifGenerationMatch'] == '0'
        self.insert_calls += 1

        def execute():
            self.manifest = self.collision or kwargs['body']
            if self.lost_write or self.collision:
                msg = 'lost response'
                raise OSError(msg)
            return self.manifest

        return Mock(execute=execute)


def test_gcs_run_binding_same_retry_changed_retry():
    client = FakeGcs()
    ensure_run_binding('gs://synthetic-output/run', {'identity': IDENTITY}, client)
    ensure_run_binding('gs://synthetic-output/run', {'identity': IDENTITY}, client)
    assert client.insert_calls == 1
    with pytest.raises(ValueError, match='different'):
        ensure_run_binding(
            'gs://synthetic-output/run',
            {'identity': IDENTITY, 'parameter': 'changed'},
            client,
        )


def test_gcs_existing_unbound_run_rejected():
    client = FakeGcs(old_outputs=True)
    with pytest.raises(ValueError, match='unbound'):
        ensure_run_binding('gs://synthetic-output/run', {'identity': IDENTITY}, client)
    assert client.insert_calls == 0


def test_gcs_lost_success_reconciled_without_rewrite():
    client = FakeGcs(lost_write=True)
    ensure_run_binding('gs://synthetic-output/run', {'identity': IDENTITY}, client)
    assert client.insert_calls == 1


def test_gcs_concurrent_different_binding_rejected():
    client = FakeGcs(collision={'metadata': {'binding_sha256': 'different'}})
    with pytest.raises(ValueError, match='different'):
        ensure_run_binding('gs://synthetic-output/run', {'identity': IDENTITY}, client)
    assert client.insert_calls == 1


def test_gcs_lost_uncommitted_write_fails_closed():
    client = FakeGcs()
    client.insert = Mock(side_effect=OSError('not committed'))
    with pytest.raises(OSError, match='not committed'):
        ensure_run_binding('gs://synthetic-output/run', {'identity': IDENTITY}, client)


def test_run_binding_scientific_parameters_and_retry():
    from types import SimpleNamespace

    from loading_pipeline.lib.misc.callset_identity import run_binding

    task = SimpleNamespace(
        callset_path=PATH,
        callset_identity=IDENTITY,
        reference_genome='GRCh38',
        dataset_type='SNV_INDEL',
        sample_type='WGS',
        project_guids=['synthetic_b', 'synthetic_a'],
        validations_to_skip=[],
        skip_check_sex_and_relatedness=False,
        skip_expect_tdr_metrics=False,
        attempt_id=0,
    )
    initial = run_binding(task)
    task.attempt_id = 1
    task.project_guids.reverse()
    assert run_binding(task) == initial
    task.skip_expect_tdr_metrics = True
    assert run_binding(task) != initial


def test_request_serialization_revalidates_identity():
    import importlib.util
    import sys
    from enum import Enum
    from pathlib import Path
    from types import ModuleType

    from pydantic import ValidationError

    core = ModuleType('loading_pipeline.lib.core')
    core.DatasetType = Enum('DatasetType', {'SNV_INDEL': 'SNV_INDEL'})
    core.ReferenceGenome = Enum('ReferenceGenome', {'GRCh38': 'GRCh38'})
    core.SampleType = Enum('SampleType', {'WGS': 'WGS'})
    clickhouse = ModuleType('loading_pipeline.lib.misc.clickhouse')
    clickhouse.ClickhouseReferenceDataset = Enum('Reference', {'TEST': 'TEST'})
    validation = ModuleType('loading_pipeline.lib.misc.validation')
    validation.ALL_VALIDATIONS = 'all'
    validation.SKIPPABLE_VALIDATIONS = [lambda: None]
    hailtop = ModuleType('hailtop')
    fs = ModuleType('hailtop.fs')
    fs.exists = Mock(return_value=True)
    hailtop.fs = fs
    modules = {m.__name__: m for m in (core, clickhouse, validation, hailtop, fs)}
    source = Path(__file__).parents[2] / 'api' / 'model.py'
    spec = importlib.util.spec_from_file_location('synthetic_pipeline_model', source)
    module = importlib.util.module_from_spec(spec)
    with (
        patch.dict(sys.modules, modules),
        patch(
            'loading_pipeline.lib.misc.callset_identity.verify_callset_identity',
        ) as verify,
    ):
        spec.loader.exec_module(module)
        request = module.LoadingPipelineRequest(
            callset_path=PATH,
            callset_identity=IDENTITY,
            project_guids=['synthetic_project'],
            sample_type='WGS',
            reference_genome='GRCh38',
            dataset_type='SNV_INDEL',
        )
        assert request.incr_attempt()
        serialized = request.model_dump_json()
        replayed = module.LoadingPipelineRequest.model_validate_json(serialized)
        assert replayed.callset_identity == IDENTITY
        assert replayed.attempt_id == 1
        verify.side_effect = ValueError('replaced between validation and enqueue')
        with pytest.raises(ValidationError, match='replaced'):
            module.LoadingPipelineRequest.model_validate_json(serialized)
