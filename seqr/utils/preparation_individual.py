"""Restricted singleton metadata adapter; only used by an explicitly enabled worker."""
from seqr.utils.preparation_ledger import PreparationConflict, canonical_digest

METADATA_FIELDS = ('birth_year', 'death_year', 'notes', 'consanguinity', 'affected_relatives', 'onset_age', 'expected_inheritance')
TOKEN_FIELDS = ('individual_id', 'family_id', 'father_id', 'mother_id', 'sex', 'affected', 'features', 'absent_features', 'candidate_genes', *METADATA_FIELDS)


def metadata_token(individual):
    value = {field: getattr(individual, field) for field in TOKEN_FIELDS}
    # Also detects edits to fields this adapter does not write.
    value['last_modified_date'] = individual.last_modified_date.isoformat() if individual.last_modified_date else None
    return canonical_digest(value)


def apply_individual_metadata(envelope, ledger, user):
    from seqr.models import Project, Family, Individual
    from seqr.views.utils.json_to_orm_utils import create_model_from_json, update_individual_from_json

    accession = envelope['accession_id']
    data = envelope['individual']
    metadata = data['metadata']
    if data['sample_id'] != accession or envelope['genome_version'] != '38':
        raise PreparationConflict('Unsupported identity or reference')
    # Existing uploader sends free text, while the model requires structured gene IDs.
    if metadata.get('candidate_genes'):
        raise PreparationConflict('Candidate genes require validated structured mapping')
    projects = Project.objects.select_for_update().filter(name=accession)
    if projects.count() > 1:
        raise PreparationConflict('Ambiguous project identity')
    project = projects.first()
    if project and (project.genome_version != envelope['genome_version'] or (ledger.project_guid and project.guid != ledger.project_guid)):
        raise PreparationConflict('Project reference or identity changed')
    if project:
        people = list(Individual.objects.select_for_update().filter(family__project=project).select_related('family'))
        if len(people) != 1 or people[0].individual_id != accession or people[0].family.family_id != accession:
            raise PreparationConflict('Existing project is not the expected singleton')
        individual = people[0]
        if not ledger.review_sequence:
            raise PreparationConflict('Existing clinical records require explicitly reviewed baseline adoption')
        if individual.active_datasets.exists() or individual.inactive_datasets.exists():
            raise PreparationConflict('Loaded individual updates require a separate validated refresh workflow')
        if metadata_token(individual) != ledger.metadata_state_token:
            raise PreparationConflict('Clinical metadata was edited outside this preparation revision')
    else:
        if ledger.review_sequence:
            raise PreparationConflict('Previously prepared project is missing')
        project = create_model_from_json(Project, {'name': accession, 'genome_version': '38'}, user)
        family = create_model_from_json(Family, {'project': project, 'family_id': accession}, user)
        individual = create_model_from_json(Individual, {'family': family, 'individual_id': accession, 'affected': data['affected']}, user)
    updates = {'sex': data['sex'], 'affected': data['affected'],
               'features': [{'id': term} for term in metadata['hpo_present']],
               'absent_features': [{'id': term} for term in metadata['hpo_absent']]}
    for field in METADATA_FIELDS:
        if field in metadata:
            updates[field] = metadata[field] or None if field == 'onset_age' else metadata[field]
    if 'candidate_genes' in metadata:
        updates['candidate_genes'] = []
    update_individual_from_json(individual, updates, user, allow_features_update=True, allow_search_field_update=True)
    individual.refresh_from_db()
    return project.guid, metadata_token(individual)


def refresh_preparation_dictionaries(user):
    from clickhouse_search.models.postgres_dicts import AffectedDict, SexDict, IndividualMetadataDict
    for dictionary in (AffectedDict, SexDict, IndividualMetadataDict):
        dictionary.reload(user)
