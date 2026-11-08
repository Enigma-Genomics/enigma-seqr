# ruff: noqa: PLR2004
# Copyright 2026 Enigma Genomics
"""Offline, fail-closed family preparation contract. Never loads or mutates seqr."""

import argparse
import json
import re
import sys
from copy import deepcopy
from pathlib import Path


class ManifestError(ValueError):
    pass


def require(condition, message):
    if not condition:
        raise ManifestError(message)


def identifier(value):
    return isinstance(value, str) and re.fullmatch(
        r'[A-Za-z0-9][A-Za-z0-9_.-]{0,127}',
        value,
    )


def prep_plan(manifest, observed_sources):  # noqa: PLR0915 - explicit contract checks
    """Validate pinned sources against separately observed headers, without IO.

    observed_sources is keyed by accession and contains bucket/path/generation,
    sample_ids, genome_version, sample_type, dataset_type. A trusted worker must
    derive it from pinned object metadata and VCF headers, never copy the manifest.
    """
    require(isinstance(manifest, dict), 'Manifest must be an object')
    require(manifest.get('schema_version') == 2, 'Requires schema_version 2')
    require(manifest.get('kind') == 'family', 'Requires family manifest')
    require(identifier(manifest.get('family_id')), 'Invalid family_id')
    require(identifier(manifest.get('joint_report_id')), 'Invalid joint_report_id')
    require(manifest.get('genome_version') == '38', 'Only GRCh38 validated')
    require(manifest.get('sample_type') in ('WES', 'WGS'), 'Invalid sample_type')
    require(
        manifest.get('dataset_type', 'SNV_INDEL') == 'SNV_INDEL',
        'Unsupported dataset',
    )
    require(manifest.get('manual_load_required') is True, 'Manual Load Data required')
    people = manifest.get('individuals')
    require(
        isinstance(people, list) and len(people) == 3,
        'Requires exactly three people',
    )
    require(all(isinstance(p, dict) for p in people), 'Invalid individuals')
    require(
        {p.get('role') for p in people} == {'proband', 'mother', 'father'},
        'Invalid roles',
    )
    for field in ('sample_id', 'accession_id'):
        values = [p.get(field) for p in people]
        require(all(identifier(v) for v in values), f'Invalid {field}')
        require(len(set(values)) == 3, f'Duplicate {field}')
    require(isinstance(observed_sources, dict), 'Missing source observations')
    require(
        set(observed_sources) == {p['accession_id'] for p in people},
        'Observation set differs',
    )
    by_role = {p['role']: p for p in people}
    sources = set()
    rows = []
    for person in people:
        sid = person['sample_id']
        require(person.get('sex') in ('M', 'F', 'U'), f'{sid}: invalid sex')
        require(person.get('affected') in ('A', 'N', 'U'), f'{sid}: invalid affected')
        for parent in ('father', 'mother'):
            expected = (
                by_role[parent]['sample_id'] if person['role'] == 'proband' else None
            )
            require(
                f'{parent}_id' in person and person[f'{parent}_id'] == expected,
                f'{sid}: explicit {parent}_id must match pedigree',
            )
        source = person.get('source')
        require(isinstance(source, dict), f'{sid}: missing source')
        for field in ('bucket', 'path', 'generation', 'vcf_sample_id'):
            require(
                isinstance(source.get(field), str) and bool(source[field].strip()),
                f'{sid}: missing source {field}',
            )
        require(source['generation'].isdigit(), f'{sid}: invalid generation')
        key = (source['bucket'], source['path'], source['generation'])
        require(key not in sources, 'Duplicate source object')
        sources.add(key)
        observed = observed_sources[person['accession_id']]
        require(isinstance(observed, dict), f'{sid}: invalid observation')
        for field in ('bucket', 'path', 'generation'):
            require(
                observed.get(field) == source[field],
                f'{sid}: source {field} mismatch',
            )
        require(
            observed.get('sample_ids') == [source['vcf_sample_id']],
            f'{sid}: header mismatch',
        )
        for field, value in (
            ('genome_version', '38'),
            ('sample_type', manifest['sample_type']),
            ('dataset_type', 'SNV_INDEL'),
        ):
            require(observed.get(field) == value, f'{sid}: observed {field} mismatch')
            require(source.get(field) == value, f'{sid}: source {field} mismatch')
        metadata = person.get('metadata', {})
        require(isinstance(metadata, dict), f'{sid}: invalid metadata')
        for field in ('hpo_present', 'hpo_absent'):
            values = metadata.get(field, [])
            require(
                isinstance(values, list) and len(values) <= 500,
                f'{sid}: invalid {field}',
            )
            require(
                all(
                    isinstance(v, str) and re.fullmatch(r'HP:\d{7}', v) for v in values
                ),
                f'{sid}: invalid HPO',
            )
            require(len(values) == len(set(values)), f'{sid}: duplicate HPO')
        require(
            not set(metadata.get('hpo_present', []))
            & set(metadata.get('hpo_absent', [])),
            f'{sid}: conflicting HPO',
        )
        require(
            metadata.get('onset_age', '') in ('', *'GEFNICJAYML'),
            f'{sid}: invalid onset',
        )
        inheritance = metadata.get('expected_inheritance', [])
        require(
            isinstance(inheritance, list)
            and all(
                isinstance(v, str) and v in list('SDLACRGXZYWFM') for v in inheritance
            ),
            f'{sid}: invalid inheritance',
        )
        rows.append(
            {
                'family_id': manifest['family_id'],
                'individual_id': sid,
                'father_id': person['father_id'],
                'mother_id': person['mother_id'],
                'sex': person['sex'],
                'affected': person['affected'],
                'features': [{'id': h} for h in metadata.get('hpo_present', [])],
                'absent_features': [{'id': h} for h in metadata.get('hpo_absent', [])],
                'onset_age': metadata.get('onset_age') or None,
                'expected_inheritance': inheritance,
                'source_metadata': deepcopy(metadata),
            },
        )
    return {
        'status': 'validated_plan_only',
        'prepared': False,
        'searchable': False,
        'manual_load_required': True,
        'family_id': manifest['family_id'],
        'individuals': rows,
        'limitations': [
            'Final single-sample VCF merge is not joint genotyping.',
            'Absent calls must remain missing, never fill with homozygous reference.',
            'No source downloads, merge, index, clinical writes or loading performed.',
            'Kinship, sex QC and multi-sample Exomiser remain unvalidated.',
        ],
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('manifest', type=Path)
    parser.add_argument('observations', type=Path)
    args = parser.parse_args()
    try:
        plan = prep_plan(
            json.loads(args.manifest.read_text()),
            json.loads(args.observations.read_text()),
        )
    except (ManifestError, ValueError, OSError) as error:
        parser.exit(2, f'Preparation blocked: {error}\n')
    sys.stdout.write(json.dumps(plan, indent=2) + '\n')


if __name__ == '__main__':
    main()
