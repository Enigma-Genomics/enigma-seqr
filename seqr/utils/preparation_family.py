"""Atomic family metadata adapter. Never enabled implicitly by singleton workers."""

import datetime

from django.conf import settings
from django.db import connection, transaction
from django.db.models import FileField

from seqr.utils.preparation_family_contract import validate_family_envelope
from seqr.utils.preparation_ledger import PreparationConflict, canonical_digest


def _snapshot(model):
    result = {}
    for field in model._meta.concrete_fields:
        value = getattr(model, field.attname)
        if isinstance(field, FileField):
            value = value.name
        result[field.attname] = (
            value.isoformat()
            if isinstance(value, (datetime.date, datetime.datetime))
            else value
        )
    return result


def family_metadata_token(project, family, people):
    return canonical_digest(
        {
            'project': _snapshot(project),
            'family': _snapshot(family),
            'individuals': [
                _snapshot(person)
                for person in sorted(people, key=lambda p: p.individual_id)
            ],
        }
    )


@transaction.atomic
def apply_family_metadata(envelope, ledger, user):
    """Caller must hold the namespaced ledger lock; this decorator cannot prove that.

    No external IO or dataset loading occurs in this callback.
    """
    from seqr.models import Project, Family, Individual
    from seqr.views.utils.json_to_orm_utils import (
        create_model_from_json,
        update_individual_from_json,
    )

    if not getattr(settings, 'PREPARATION_FAMILY_ENABLED', False):
        raise PreparationConflict('Family preparation activation is disabled')
    subject = validate_family_envelope(envelope)
    if not connection.in_atomic_block or ledger.subject_key != subject:
        raise PreparationConflict(
            'Family adapter requires its locked namespaced ledger'
        )
    if ledger.review_sequence and ledger.request_id != envelope['request_id']:
        raise PreparationConflict('Family membership, roles, source or assay changed')
    expected_id = envelope['family_id']
    projects = list(Project.objects.select_for_update().filter(name=expected_id))
    if len(projects) > 1:
        raise PreparationConflict('Ambiguous family project identity')
    if projects:
        project = projects[0]
        if (
            not ledger.review_sequence
            or not ledger.project_guid
            or ledger.project_guid != project.guid
            or project.genome_version != '38'
        ):
            raise PreparationConflict('Cannot adopt existing or changed family project')
        families = list(Family.objects.select_for_update().filter(project=project))
        if len(families) != 1 or families[0].family_id != expected_id:
            raise PreparationConflict('Family membership changed')
        family = families[0]
        people = list(Individual.objects.select_for_update().filter(family=family))
        if len(people) != 3 or {p.individual_id for p in people} != {
            p['sample_id'] for p in envelope['individuals']
        }:
            raise PreparationConflict('Individual membership changed')
        if any(
            p.active_datasets.exists() or p.inactive_datasets.exists() for p in people
        ):
            raise PreparationConflict(
                'Loaded family requires a separate refresh workflow'
            )
        if (
            family_metadata_token(project, family, people)
            != ledger.metadata_state_token
        ):
            raise PreparationConflict(
                'Family metadata changed outside reviewed preparation'
            )
    else:
        if ledger.review_sequence or ledger.project_guid:
            raise PreparationConflict('Previously owned family project missing')
        project = create_model_from_json(
            Project, {'name': expected_id, 'genome_version': '38'}, user
        )
        family = create_model_from_json(
            Family, {'project': project, 'family_id': expected_id}, user
        )
        people = [
            create_model_from_json(
                Individual,
                {
                    'family': family,
                    'individual_id': p['sample_id'],
                    'sex': p['sex'],
                    'affected': p['affected'],
                },
                user,
            )
            for p in envelope['individuals']
        ]
    by_id = {p.individual_id: p for p in people}
    for data in envelope['individuals']:
        person = by_id[data['sample_id']]
        updates = {
            'sex': data['sex'],
            'affected': data['affected'],
            'mother': by_id.get(data['mother_id']),
            'father': by_id.get(data['father_id']),
            'features': [{'id': term} for term in data['metadata']['hpo_present']],
            'absent_features': [
                {'id': term} for term in data['metadata']['hpo_absent']
            ],
        }
        update_individual_from_json(
            person,
            updates,
            user,
            allow_features_update=True,
            allow_search_field_update=True,
        )
        person.refresh_from_db()
    project.refresh_from_db()
    family.refresh_from_db()
    return project.guid, family_metadata_token(project, family, people)
