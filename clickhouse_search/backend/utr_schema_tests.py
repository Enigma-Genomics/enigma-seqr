from types import SimpleNamespace
from unittest.mock import Mock, patch

from django.test import SimpleTestCase

from clickhouse_search.backend import utr_schema


class UtrSchemaSafetyTests(SimpleTestCase):
    def _editor(self, responses):
        editor = Mock()
        editor.connection.cursor.return_value.fetchone.side_effect = responses
        return editor

    def test_unknown_schema_addition_is_rejected(self):
        old = SimpleNamespace(base_fields=[])
        new = SimpleNamespace(base_fields=[('unexpected', Mock())])
        with self.assertRaisesMessage(ValueError, 'Unexpected UTR schema addition: unexpected'):
            utr_schema._tuple_expression('row', old, new, Mock())

    def test_non_atomic_database_is_rejected_before_changes(self):
        editor = self._editor([('Ordinary',)])
        with self.assertRaisesMessage(RuntimeError, 'requires an Atomic database'):
            utr_schema.migrate_utr_tables(Mock(get_model=Mock(return_value=SimpleNamespace(_meta=SimpleNamespace(db_table='variants')))), editor)
        editor.create_model.assert_not_called()
        self.assertEqual(editor.connection.cursor.return_value.execute.call_count, 1)

    def test_populated_database_requires_explicit_maintenance_gate(self):
        editor = self._editor([('Atomic',), (1,)])
        with patch.dict(utr_schema.os.environ, {}, clear=True):
            with self.assertRaisesMessage(RuntimeError, 'requires verified backups and quiescent writers'):
                utr_schema.migrate_utr_tables(Mock(get_model=Mock(return_value=SimpleNamespace(_meta=SimpleNamespace(db_table='variants')))), editor)
        editor.create_model.assert_not_called()

    def test_existing_rollback_table_is_never_overwritten(self):
        editor = self._editor([('Atomic',), (0,), (1,)])
        with self.assertRaisesMessage(RuntimeError, 'stage/rollback table already exists'):
            utr_schema.migrate_utr_tables(Mock(get_model=Mock(return_value=SimpleNamespace(_meta=SimpleNamespace(db_table='variants')))), editor)
        editor.create_model.assert_not_called()

    def _migration(self, signatures, fail_second_exchange=False):
        editor = self._editor([('Atomic',), (0,), (0,), (0,), (0,), *signatures])
        field = SimpleNamespace(name='transcripts', column='transcripts')
        target = SimpleNamespace(_meta=SimpleNamespace(local_fields=[field], get_field=lambda _: field))
        original = SimpleNamespace(_meta=SimpleNamespace(
            db_table='variants', local_fields=[field], get_field=lambda _: field,
        ))
        apps = Mock()
        apps.get_model.side_effect = [original, SimpleNamespace(_meta=SimpleNamespace(
            db_table='details', local_fields=[field], get_field=lambda _: field,
        ))]
        state = Mock(name='state')
        state.name = 'Stage'
        state.options = {}
        state.fields = {}
        state.render.return_value = target
        cursor = editor.connection.cursor.return_value
        if fail_second_exchange:
            def execute(sql):
                if sql == 'EXCHANGE TABLES `details` AND `details__before_utr_v1`':
                    raise RuntimeError('simulated exchange failure')
            cursor.execute.side_effect = execute
        with patch.object(utr_schema, 'TABLES', [('One', 'transcripts', '/synthetic/one'), ('Two', 'transcripts', '/synthetic/two')]), \
                patch.object(utr_schema.ModelState, 'from_model', return_value=state), \
                patch.object(utr_schema, '_expanded_field', return_value=field), \
                patch.object(utr_schema, '_tuple_expression', return_value='converted'):
            with self.assertRaises(RuntimeError):
                utr_schema.migrate_utr_tables(apps, editor)
        return [call.args[0] for call in cursor.execute.call_args_list]

    def test_copy_signature_mismatch_preserves_original_and_drops_stage(self):
        statements = self._migration([(1, 1, 100), (1, 1, 101)])
        self.assertFalse(any(sql.startswith('EXCHANGE') for sql in statements))
        self.assertEqual(statements[-1], 'DROP TABLE `variants__before_utr_v1`')

    def test_partial_exchange_failure_rolls_back_before_stage_cleanup(self):
        statements = self._migration([(1, 1, 100)] * 4, fail_second_exchange=True)
        self.assertEqual(statements[-3:], [
            'EXCHANGE TABLES `variants` AND `variants__before_utr_v1`',
            'DROP TABLE `variants__before_utr_v1`',
            'DROP TABLE `details__before_utr_v1`',
        ])
        self.assertEqual(statements.count('EXCHANGE TABLES `variants` AND `variants__before_utr_v1`'), 2)
