# SQL fixtures use only hard-coded synthetic values and test table names.
# ruff: noqa: S608

from django.db import connections

from loading_pipeline.lib.core.environment import Env
from loading_pipeline.lib.misc.key_integrity import (
    VariantKeyIntegrityError,
    validate_entry_variant_keys,
    validate_new_variant_keys,
)
from loading_pipeline.lib.test.clickhouse_schema_testcase import (
    ClickhouseSchemaTestCase,
)


class VariantKeyIntegrityTest(ClickhouseSchemaTestCase):
    def setUp(self):
        super().setUp()
        for name, columns in [
            ('source_variants', 'key Nullable(Int64), variantId Nullable(String)'),
            ('lookup', 'key UInt32, variantId String'),
            ('source_entries', 'variantId Nullable(String)'),
        ]:
            with connections['clickhouse_write'].cursor() as cursor:
                cursor.execute(
                    f'CREATE TABLE {self.table(name)} ({columns}) ENGINE=Memory',
                )

    def tearDown(self):
        for name in ['source_variants', 'lookup', 'source_entries']:
            with connections['clickhouse_write'].cursor() as cursor:
                cursor.execute(f'DROP TABLE IF EXISTS {self.table(name)}')
        super().tearDown()

    @staticmethod
    def table(name):
        return f'{Env.CLICKHOUSE_DATABASE}.`key_integrity_{name}`'

    def execute(self, sql):
        self.assertTrue(sql.lstrip().startswith(('SELECT', 'WITH')))
        with connections['clickhouse_write'].cursor() as cursor:
            cursor.execute(sql)
            return cursor.fetchall()

    def insert(self, table, values):
        with connections['clickhouse_write'].cursor() as cursor:
            cursor.execute(f'INSERT INTO {self.table(table)} VALUES {values}')

    def validate_variants(self):
        validate_new_variant_keys(
            self.execute, self.table('source_variants'), self.table('lookup'),
        )

    def validate_entries(self):
        validate_entry_variant_keys(
            self.execute, self.table('source_entries'), self.table('lookup'),
        )

    def test_new_mapping(self):
        self.insert('lookup', "(1,'A')")
        self.insert('source_variants', "(2,'B')")
        self.validate_variants()

    def test_identical_reload(self):
        self.insert('lookup', "(1,'A')")
        self.insert('source_variants', "(1,'A')")
        self.validate_variants()

    def test_existing_key_collision(self):
        self.insert('lookup', "(1,'A')")
        self.insert('source_variants', "(1,'B')")
        with self.assertRaisesRegex(VariantKeyIntegrityError, 'different variant'):
            self.validate_variants()

    def test_existing_variant_moves_key(self):
        self.insert('lookup', "(1,'A')")
        self.insert('source_variants', "(2,'A')")
        with self.assertRaisesRegex(VariantKeyIntegrityError, 'different key'):
            self.validate_variants()

    def test_duplicate_source_key(self):
        self.insert('source_variants', "(1,'A'),(1,'B')")
        with self.assertRaisesRegex(VariantKeyIntegrityError, 'not one-to-one'):
            self.validate_variants()

    def test_duplicate_source_variant(self):
        self.insert('source_variants', "(1,'A'),(2,'A')")
        with self.assertRaisesRegex(VariantKeyIntegrityError, 'not one-to-one'):
            self.validate_variants()

    def test_invalid_source_key(self):
        self.insert('source_variants', "(-1,'A')")
        with self.assertRaisesRegex(VariantKeyIntegrityError, 'invalid source'):
            self.validate_variants()

    def test_null_source_identifier(self):
        self.insert('source_variants', '(1,NULL)')
        with self.assertRaisesRegex(VariantKeyIntegrityError, 'invalid source'):
            self.validate_variants()

    def test_uint32_overflow(self):
        self.insert('source_variants', "(4294967296,'A')")
        with self.assertRaisesRegex(VariantKeyIntegrityError, 'invalid source'):
            self.validate_variants()

    def test_empty_source(self):
        self.validate_variants()
        self.validate_entries()

    def test_valid_entries(self):
        self.insert('lookup', "(1,'A'),(2,'B')")
        self.insert('source_entries', "('A'),('A'),('B')")
        self.validate_entries()

    def test_missing_entry_mapping(self):
        self.insert('source_entries', "('A')")
        with self.assertRaisesRegex(VariantKeyIntegrityError, 'exactly one key'):
            self.validate_entries()

    def test_ambiguous_entry_key(self):
        self.insert('lookup', "(1,'A'),(1,'B')")
        self.insert('source_entries', "('A')")
        with self.assertRaisesRegex(VariantKeyIntegrityError, 'multiple variants'):
            self.validate_entries()

    def test_entry_variant_multiple_keys(self):
        self.insert('lookup', "(1,'A'),(2,'A')")
        self.insert('source_entries', "('A')")
        with self.assertRaisesRegex(VariantKeyIntegrityError, 'exactly one key'):
            self.validate_entries()

    def test_unrelated_collision(self):
        self.insert('lookup', "(1,'A'),(2,'B'),(2,'C')")
        self.insert('source_entries', "('A')")
        self.validate_entries()

    def test_null_entry(self):
        self.insert('source_entries', '(NULL)')
        with self.assertRaisesRegex(VariantKeyIntegrityError, 'no variant identifier'):
            self.validate_entries()
