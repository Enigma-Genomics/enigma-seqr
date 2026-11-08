"""Synthetic-only family preparation envelope fixture."""

from seqr.utils.preparation_family_contract import (
    family_id,
    family_subject_key,
    family_request,
)
from seqr.utils.preparation_ledger import canonical_digest


def seal(envelope):
    envelope['revision'] = canonical_digest(
        {k: v for k, v in envelope.items() if k != 'revision'}
    )
    return envelope


def synthetic_family_envelope():
    joint = 'SYNTHETIC_JOINT_REPORT'
    people = []
    for role in ('father', 'proband', 'mother'):
        people.append(
            dict(
                role=role,
                sample_id='SYNTH_' + role,
                accession_id='ACC_' + role,
                source_sample_id='VCF_' + role,
                mother_id='SYNTH_mother' if role == 'proband' else None,
                father_id='SYNTH_father' if role == 'proband' else None,
                sex='U',
                affected='U',
                metadata=dict(hpo_present=['HP:0000001'], hpo_absent=[]),
            )
        )
    envelope = dict(
        schema_version=3,
        kind='family_joint_callset',
        joint_report_id=joint,
        family_id=family_id(joint),
        subject_key=family_subject_key(joint),
        genome_version='38',
        sample_type='WES',
        dataset_type='SNV_INDEL',
        source=dict(bucket='synthetic', path='fixture.vcf.gz', generation='1'),
        source_sha256='a' * 64,
        vcf_sha256='b' * 64,
        index_sha256='c' * 64,
        records=2,
        joint_genotyping_attestation=dict(
            confirmed=True,
            authority='user_confirmed',
            reference='synthetic-run-record',
            confirmed_by='synthetic-operator',
        ),
        metadata_review=dict(
            confirmed=True,
            reference='synthetic-review-1',
            reviewer='synthetic-reviewer',
        ),
        calling_history_independently_verified=False,
        clinical_ready=False,
        prepared=False,
        searchable=False,
        manual_load_required=True,
        status='awaiting_consumer',
        review_sequence=1,
        expected_previous_revision=None,
        individuals=people,
    )
    envelope['normalized_outputs'] = {
        kind: dict(
            bucket='synthetic',
            path='normalized.' + extension,
            generation='1',
            sha256=digest * 64,
        )
        for kind, extension, digest in (
            ('vcf', 'vcf.gz', 'b'),
            ('index', 'vcf.gz.tbi', 'c'),
        )
    }
    envelope['request_id'] = canonical_digest(family_request(envelope))
    return seal(envelope)
