import json
from pathlib import Path

from django.apps import apps
from django.test import SimpleTestCase


class ClickhouseTranscriptFixtureSchemaTests(SimpleTestCase):
    def test_transcript_fixture_values_match_model_fields(self):
        fixture_dir = Path(__file__).parent / 'fixtures'
        fixture_paths = sorted(fixture_dir.glob('*.json')) + sorted(
            (fixture_dir.parent.parent / 'loading_pipeline' / 'lib' / 'test' / 'fixtures').glob('*.json')
        )
        fixtures = [row for path in fixture_paths for row in json.loads(path.read_text())]
        transcript_rows = 0
        for row in fixtures:
            field_name = 'transcripts' if 'transcripts' in row['fields'] else 'sorted_transcript_consequences'
            values = row['fields'].get(field_name)
            if values is None:
                continue
            field = apps.get_model(row['model'])._meta.get_field(field_name)
            for value in values:
                with self.subTest(model=row['model'], pk=row['pk'], value=value):
                    self.assertEqual(len(value), len(field.base_fields))
                    # Exercise the same conversion used by Django's fixture loader.
                    self.assertEqual(len(field.to_python([value])), 1)
                    transcript_rows += 1
        self.assertGreater(transcript_rows, 0)
