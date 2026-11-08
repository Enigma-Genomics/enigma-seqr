import json
from pathlib import Path

from django.apps import apps
from django.test import SimpleTestCase


class ClickhouseTranscriptFixtureSchemaTests(SimpleTestCase):
    def test_transcript_fixture_values_match_model_fields(self):
        fixture_path = Path(__file__).parent / 'fixtures' / 'clickhouse_search.json'
        fixtures = json.loads(fixture_path.read_text())
        transcript_rows = 0
        for row in fixtures:
            values = row['fields'].get('sorted_transcript_consequences')
            if values is None:
                continue
            field = apps.get_model(row['model'])._meta.get_field('sorted_transcript_consequences')
            for value in values:
                with self.subTest(model=row['model'], pk=row['pk'], value=value):
                    self.assertEqual(len(value), len(field.base_fields))
                    # Exercise the same conversion used by Django's fixture loader.
                    self.assertEqual(len(field.to_python([value])), 1)
                    transcript_rows += 1
        self.assertGreater(transcript_rows, 0)
