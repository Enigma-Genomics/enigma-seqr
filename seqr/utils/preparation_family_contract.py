"""Separate, disabled-by-default joint-callset family preparation envelope."""

import hashlib
import re

from seqr.utils.preparation_ledger import PreparationConflict, canonical_digest

ROLES = ('proband', 'mother', 'father')
IDENTITY_FIELDS = (
    'role',
    'sample_id',
    'accession_id',
    'source_sample_id',
    'mother_id',
    'father_id',
)


def require(condition, message):
    if not condition:
        raise PreparationConflict(message)


def identifier(value):
    return isinstance(value, str) and re.fullmatch(
        r'[A-Za-z0-9][A-Za-z0-9_.-]{0,127}', value
    )


def identity_hash(joint_report_id):
    require(identifier(joint_report_id), 'Invalid authoritative joint report identity')
    return hashlib.sha256(joint_report_id.encode('utf-8')).hexdigest()


def family_id(joint_report_id):
    return 'JR_' + identity_hash(joint_report_id)


def family_subject_key(joint_report_id):
    return 'family:' + identity_hash(joint_report_id)


def family_request(envelope):
    """Canonical digest of this object freezes all non-metadata inputs across reviews."""
    people = {p['role']: p for p in envelope['individuals']}
    return {
        'schema_version': 1,
        'kind': 'normalize_joint_trio',
        'joint_report_id': envelope['joint_report_id'],
        'family_id': envelope['family_id'],
        'intake_request_id': envelope['intake_request_id'],
        'genome_version': envelope['genome_version'],
        'sample_type': envelope['sample_type'],
        'dataset_type': envelope['dataset_type'],
        'source': envelope['source'],
        'source_sha256': envelope['source_sha256'],
        'normalized_outputs': envelope['normalized_outputs'],
        'joint_genotyping_attestation': envelope['joint_genotyping_attestation'],
        'individuals': [
            {key: people[role][key] for key in IDENTITY_FIELDS} for role in ROLES
        ],
    }


def validate_family_envelope(envelope):
    require(isinstance(envelope, dict), 'Envelope must be an object')
    require(
        envelope.get('schema_version') == 3
        and envelope.get('kind') == 'family_joint_callset',
        'Unsupported family schema',
    )
    require(
        envelope.get('revision')
        == canonical_digest({k: v for k, v in envelope.items() if k != 'revision'}),
        'Revision digest mismatch',
    )
    require(
        envelope.get('family_id') == family_id(envelope.get('joint_report_id')),
        'Derived family identity mismatch',
    )
    require(
        envelope.get('subject_key') == family_subject_key(envelope['joint_report_id']),
        'Family subject namespace mismatch',
    )
    require(
        envelope.get('genome_version') == '38'
        and envelope.get('sample_type') in ('WES', 'WGS')
        and envelope.get('dataset_type') == 'SNV_INDEL',
        'Unsupported assay/build',
    )
    require(
        envelope.get('status') == 'awaiting_consumer'
        and envelope.get('prepared') is False
        and envelope.get('searchable') is False
        and envelope.get('manual_load_required') is True,
        'Invalid lifecycle',
    )
    require(
        envelope.get('clinical_ready') is False
        and envelope.get('calling_history_independently_verified') is False,
        'Unsupported clinical/provenance certification',
    )
    sequence, previous = (
        envelope.get('review_sequence'),
        envelope.get('expected_previous_revision'),
    )
    require(
        type(sequence) is int and 1 <= sequence <= 9007199254740991,
        'Invalid review sequence',
    )
    require(
        (sequence == 1 and previous is None)
        or (
            sequence > 1
            and isinstance(previous, str)
            and re.fullmatch('[a-f0-9]{64}', previous)
        ),
        'Invalid predecessor',
    )
    for key in ('source_sha256', 'vcf_sha256', 'index_sha256', 'request_id', 'intake_request_id'):
        require(
            isinstance(envelope.get(key), str)
            and re.fullmatch('[a-f0-9]{64}', envelope[key]),
            f'Invalid {key}',
        )
    require(
        type(envelope.get('records')) is int and envelope['records'] >= 0,
        'Invalid record count',
    )
    source = envelope.get('source')
    require(
        isinstance(source, dict) and set(source) == {'bucket', 'path', 'generation'},
        'Invalid pinned joint source',
    )
    for key, bound in (('bucket', 222), ('path', 1024), ('generation', 30)):
        require(
            isinstance(source[key], str)
            and 0 < len(source[key]) <= bound
            and not any(ord(c) < 32 for c in source[key]),
            'Invalid source pin',
        )
    require(
        re.fullmatch('[1-9][0-9]*', source['generation']), 'Invalid source generation'
    )
    outputs = envelope.get('normalized_outputs')
    require(
        isinstance(outputs, dict) and set(outputs) == {'vcf', 'index'},
        'Pinned normalized outputs required',
    )
    for kind, hash_key in (('vcf', 'vcf_sha256'), ('index', 'index_sha256')):
        output = outputs[kind]
        require(
            isinstance(output, dict)
            and set(output) == {'bucket', 'path', 'generation', 'sha256'},
            'Invalid normalized output pin',
        )
        for key, bound in (('bucket', 222), ('path', 1024), ('generation', 30)):
            require(
                isinstance(output[key], str)
                and 0 < len(output[key]) <= bound
                and not any(ord(c) < 32 for c in output[key]),
                'Invalid output location',
            )
        require(
            re.fullmatch('[1-9][0-9]*', output['generation'])
            and output['sha256'] == envelope[hash_key],
            'Invalid output generation/hash',
        )
    require(
        (outputs['vcf']['bucket'], outputs['vcf']['path'])
        != (outputs['index']['bucket'], outputs['index']['path']),
        'VCF/index objects must differ',
    )
    provenance = envelope.get('joint_genotyping_attestation')
    require(
        isinstance(provenance, dict)
        and set(provenance) == {'confirmed', 'authority', 'reference', 'confirmed_by'},
        'Invalid attestation',
    )
    require(
        provenance['confirmed'] is True
        and provenance['authority'] in ('user_confirmed', 'trusted_workflow'),
        'Joint provenance confirmation required',
    )
    for key in ('reference', 'confirmed_by'):
        require(
            isinstance(provenance[key], str)
            and 0 < len(provenance[key]) <= 500
            and provenance[key].strip() == provenance[key],
            'Invalid attestation reference',
        )
    people = envelope.get('individuals')
    require(
        isinstance(people, list)
        and len(people) == 3
        and all(isinstance(p, dict) for p in people),
        'Exactly three people required',
    )
    require(
        all(p.get('role') in ROLES for p in people)
        and len({p['role'] for p in people}) == 3,
        'Explicit unique roles required',
    )
    for key in ('sample_id', 'accession_id', 'source_sample_id'):
        require(
            all(identifier(p.get(key)) for p in people)
            and len({p[key] for p in people}) == 3,
            f'Invalid {key}',
        )
    by_role = {p['role']: p for p in people}
    for person in people:
        for parent in ('mother', 'father'):
            require(
                parent + '_id' in person
                and person[parent + '_id']
                == (
                    by_role[parent]['sample_id']
                    if person['role'] == 'proband'
                    else None
                ),
                'Explicit parent mapping mismatch',
            )
        require(
            person.get('sex') in ('M', 'F', 'U')
            and person.get('affected') in ('A', 'N', 'U'),
            'Invalid reviewed sex/affected',
        )
        metadata = person.get('metadata')
        require(
            isinstance(metadata, dict)
            and set(metadata) == {'hpo_present', 'hpo_absent'},
            'Only explicitly reviewed HPO metadata is supported',
        )
        for key in ('hpo_present', 'hpo_absent'):
            values = metadata[key]
            require(
                isinstance(values, list)
                and len(values) <= 500
                and all(
                    isinstance(v, str) and re.fullmatch(r'HP:\d{7}', v) for v in values
                ),
                'Invalid HPO list',
            )
            require(len(set(values)) == len(values), 'Duplicate HPO term')
        require(
            not set(metadata['hpo_present']) & set(metadata['hpo_absent']),
            'Conflicting HPO terms',
        )
    review = envelope.get('metadata_review')
    require(
        isinstance(review, dict)
        and set(review) == {'confirmed', 'reference', 'reviewer'}
        and review['confirmed'] is True,
        'Separate metadata review required',
    )
    for key in ('reference', 'reviewer'):
        require(
            isinstance(review[key], str)
            and 0 < len(review[key]) <= 500
            and bool(review[key].strip()),
            'Invalid metadata reviewer/reference',
        )
    require(
        envelope['request_id'] == canonical_digest(family_request(envelope)),
        'Family request fingerprint mismatch',
    )
    return envelope['subject_key']
