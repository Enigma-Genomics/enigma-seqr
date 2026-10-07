"""Copy and verify UTR schema changes before swapping any live table.

EmbeddedRocksDB values require a table rebuild when nested tuple schemas change.
Keep the pre-migration tables for rollback. This operation requires quiescent writes.
"""
import os
from django.db.migrations.state import ModelState
from clickhouse_backend import models
from clickhouse_search.backend.engines import EmbeddedRocksDB
from settings import CLICKHOUSE_DATA_DIR, CLICKHOUSE_IN_MEMORY_DIR

TABLES = [
    ('VariantsSnvIndel', 'sorted_transcript_consequences', CLICKHOUSE_IN_MEMORY_DIR + '/GRCh38/SNV_INDEL/variants_utr_v1'),
    ('VariantsDiskSnvIndel', 'sorted_transcript_consequences', CLICKHOUSE_DATA_DIR + '/GRCh38/SNV_INDEL/variants_utr_v1'),
    ('VariantDetailsSnvIndel', 'transcripts', CLICKHOUSE_DATA_DIR + '/GRCh38/SNV_INDEL/variants_details_utr_v1'),
]


def _expanded_field(field, details=False):
    _, _, args, kwargs = field.deconstruct()
    base = list(kwargs['base_fields'])
    if details:
        for i, (name, child) in enumerate(base):
            if name == 'utrannotator':
                _, _, child_args, child_kwargs = child.deconstruct()
                child_kwargs['base_fields'] = list(child_kwargs['base_fields']) + [
                    ('fiveutrConsequences', models.ArrayField(models.StringField())),
                    ('fiveutrEffectsJson', models.StringField(null=True, blank=True)),
                ]
                base[i] = (name, type(child)(*child_args, **child_kwargs))
    else:
        base += [('fiveutrConsequences', models.ArrayField(models.StringField()))]
        base.sort(key=lambda item: item[0])
    kwargs['base_fields'] = base
    return type(field)(*args, **kwargs)


def _quoted(value):
    return '`' + value.replace('`', '``') + '`'


def _type_literal(value):
    return "'" + value.replace('\\', '\\\\').replace("'", "\\'") + "'"


def _tuple_expression(source, old, new, connection):
    old_fields = dict(old.base_fields)
    expressions = []
    for name, field in new.base_fields:
        value = source + '.' + _quoted(name)
        if name not in old_fields:
            if name == 'fiveutrConsequences':
                previous = source + '.fiveutrConsequence'
                value = f'if(isNull({previous}), CAST([], \'Array(String)\'), [toString({previous})])'
            elif name == 'fiveutrEffectsJson':
                value = "CAST(NULL, 'Nullable(String)')"
            else:
                raise ValueError('Unexpected UTR schema addition: ' + name)
        elif hasattr(field, 'base_fields'):
            value = _tuple_expression(value, old_fields[name], field, connection)
        expressions.append(value)
    kind = new.db_type(connection)
    if kind.startswith('Nested('):
        kind = 'Tuple(' + kind[len('Nested('):]
    return f'CAST(tuple({", ".join(expressions)}), {_type_literal(kind)})'


def migrate_utr_tables(apps, schema_editor):
    connection = schema_editor.connection
    cursor = connection.cursor()
    cursor.execute('SELECT engine FROM system.databases WHERE name = currentDatabase()')
    if cursor.fetchone()[0] != 'Atomic':
        raise RuntimeError('UTR migration requires an Atomic database for table exchange')
    plans = []
    for model_name, column_name, path in TABLES:
        original = apps.get_model('clickhouse_search', model_name)
        table = original._meta.db_table
        cursor.execute(f'SELECT count() FROM {_quoted(table)}')
        if cursor.fetchone()[0] and os.environ.get('ENIGMA_UTR_SCHEMA_MIGRATION') != 'allow':
            raise RuntimeError('UTR migration requires verified backups and quiescent writers; set ENIGMA_UTR_SCHEMA_MIGRATION=allow only for the reviewed maintenance step')
        stage = table + '__before_utr_v1'
        cursor.execute('EXISTS TABLE ' + _quoted(stage))
        if cursor.fetchone()[0]:
            raise RuntimeError('UTR migration stage/rollback table already exists: ' + stage)
        state = ModelState.from_model(original)
        state.name += 'UtrStage'
        state.options['db_table'] = stage
        state.options['engine'] = EmbeddedRocksDB(0, path, primary_key='key', flatten_nested=0)
        old_field = original._meta.get_field(column_name)
        state.fields[column_name] = _expanded_field(old_field, model_name == 'VariantDetailsSnvIndel')
        target = state.render(apps)
        new_field = target._meta.get_field(column_name)
        selected = []
        for field in original._meta.local_fields:
            value = _quoted(field.column)
            if field.name == column_name:
                value = 'arrayMap(utr_row -> ' + _tuple_expression('utr_row', old_field, new_field, connection) + ', ' + value + ')'
            selected.append(value)
        plans.append((table, stage, target, selected))
    created, exchanged = [], []
    try:
        # Prepare and verify all tables before the first swap. Originals stay intact.
        for table, stage, target, selected in plans:
            schema_editor.create_model(target)
            created.append(stage)
            columns = ', '.join(_quoted(f.column) for f in target._meta.local_fields)
            cursor.execute(f'INSERT INTO {_quoted(stage)} ({columns}) SELECT {", ".join(selected)} FROM {_quoted(table)}')
            expected = f'SELECT count(), uniqExact(key), sum(cityHash64(toJSONString(tuple({", ".join(selected)})))) FROM {_quoted(table)}'
            cursor.execute(expected)
            source_signature = cursor.fetchone()
            cursor.execute(f'SELECT count(), uniqExact(key), sum(cityHash64(toJSONString(tuple({columns})))) FROM {_quoted(stage)}')
            if cursor.fetchone() != source_signature:
                raise RuntimeError('UTR table copy verification failed; original retained: ' + table)
        for table, stage, _, _ in plans:
            cursor.execute(f'EXCHANGE TABLES {_quoted(table)} AND {_quoted(stage)}')
            exchanged.append((table, stage))
    except Exception:
        for table, stage in reversed(exchanged):
            cursor.execute(f'EXCHANGE TABLES {_quoted(table)} AND {_quoted(stage)}')
        for stage in created:
            cursor.execute('DROP TABLE ' + _quoted(stage))
        raise
    # Originals remain in __before_utr_v1. Do not delete them automatically.
