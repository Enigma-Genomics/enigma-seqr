"""Transactional revision CAS and replayable acknowledgment, with no cloud IO."""
import hashlib
import json
import re

from django.db import transaction
from django.utils import timezone

from seqr.preparation_models import PreparationLedger, PreparationOutbox


class PreparationConflict(ValueError):
    pass


def canonical_digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':'), ensure_ascii=False).encode()).hexdigest()


def validate_envelope(envelope):
    revision = envelope.get('revision')
    content = {key: value for key, value in envelope.items() if key != 'revision'}
    if not isinstance(revision, str) or not re.fullmatch('[a-f0-9]{64}', revision) or canonical_digest(content) != revision:
        raise PreparationConflict('Revision does not match canonical envelope')
    if envelope.get('schema_version') != 2 or envelope.get('kind') != 'individual':
        raise PreparationConflict('Unsupported preparation schema')
    if not isinstance(envelope.get('accession_id'), str) or not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]{0,127}', envelope['accession_id']):
        raise PreparationConflict('Invalid subject identity')
    if not isinstance(envelope.get('request_id'), str) or not re.fullmatch('[a-f0-9]{64}', envelope['request_id']):
        raise PreparationConflict('Invalid request identity')
    source = envelope.get('source')
    if not isinstance(source, dict) or not all(isinstance(source.get(k), str) and source[k] for k in ('bucket', 'path', 'generation', 'vcf_sample_id')) or not re.fullmatch(r'[1-9]\d*', source['generation']):
        raise PreparationConflict('Invalid pinned source')
    request = {key: envelope[key] for key in ('accession_id', 'sample_type', 'genome_version')}
    request = {'schema_version': 1, 'kind': 'normalize_singleton', **request, 'source': {key: source[key] for key in ('bucket', 'path', 'generation', 'vcf_sample_id')}}
    request_id = hashlib.sha256(json.dumps(request, separators=(',', ':'), ensure_ascii=False).encode()).hexdigest()
    if request_id != envelope['request_id']:
        raise PreparationConflict('Source request fingerprint mismatch')
    if envelope.get('status') != 'awaiting_consumer' or type(envelope.get('records')) is not int or envelope['records'] < 0 or not isinstance(envelope.get('source_sha256'), str) or not re.fullmatch('[a-f0-9]{64}', envelope['source_sha256']):
        raise PreparationConflict('Invalid normalization evidence')
    sequence = envelope.get('review_sequence')
    previous = envelope.get('expected_previous_revision')
    if type(sequence) is not int or not 1 <= sequence <= 9007199254740991:
        raise PreparationConflict('Invalid ordered review sequence')
    if (sequence == 1 and previous is not None) or (sequence > 1 and (
            not isinstance(previous, str) or not re.fullmatch('[a-f0-9]{64}', previous))):
        raise PreparationConflict('Invalid previous revision')
    if envelope.get('genome_version') != '38' or envelope.get('sample_type') not in ('WES', 'WGS'):
        raise PreparationConflict('Unsupported dataset provenance')
    individual = envelope.get('individual', {})
    if not isinstance(individual, dict) or individual.get('sample_id') != envelope['accession_id'] or individual.get('sex') not in ('M', 'F', 'U') or individual.get('affected') not in ('A', 'N', 'U'):
        raise PreparationConflict('Invalid per-person identity or status')
    metadata = individual.get('metadata')
    if not isinstance(metadata, dict) or set(metadata) - {'hpo_present', 'hpo_absent', 'birth_year', 'death_year', 'notes', 'consanguinity', 'affected_relatives', 'candidate_genes', 'onset_age', 'expected_inheritance'}:
        raise PreparationConflict('Unsupported individual metadata')
    for key in ('hpo_present', 'hpo_absent'):
        values = metadata.get(key)
        if not isinstance(values, list) or len(values) > 500 or not all(isinstance(v, str) and re.fullmatch(r'HP:\d{7}', v) for v in values) or len(set(values)) != len(values):
            raise PreparationConflict('Invalid HPO terms')
    if set(metadata['hpo_present']) & set(metadata['hpo_absent']):
        raise PreparationConflict('Conflicting HPO terms')
    for key in ('birth_year', 'death_year'):
        if key in metadata and (type(metadata[key]) is not int or not 1900 <= metadata[key] <= 2100):
            raise PreparationConflict('Invalid year')
    for key in ('consanguinity', 'affected_relatives'):
        if key in metadata and type(metadata[key]) is not bool:
            raise PreparationConflict('Invalid boolean metadata')
    for key, limit in (('notes', 2000), ('candidate_genes', 500)):
        if key in metadata and (not isinstance(metadata[key], str) or len(metadata[key]) > limit):
            raise PreparationConflict('Invalid text metadata')
    if metadata.get('onset_age', '') not in ('', *'GEFNICJAYML'):
        raise PreparationConflict('Invalid onset')
    inheritance = metadata.get('expected_inheritance', [])
    if not isinstance(inheritance, list) or not all(isinstance(v, str) and v in list('SDLACRGXZYWFM') for v in inheritance):
        raise PreparationConflict('Invalid inheritance')
    if envelope.get('prepared') is not False or envelope.get('searchable') is not False or envelope.get('manual_load_required') is not True:
        raise PreparationConflict('Invalid preparation lifecycle')


@transaction.atomic
def apply_preparation(envelope, apply_metadata):
    """apply_metadata runs inside this transaction and returns (project_guid, state_token).

    It MUST lock existing clinical rows and compare their current metadata token to
    ledger.metadata_state_token before updates. Initial adoption of existing rows is
    prohibited without a separately reviewed baseline. It must not do external IO.
    """
    validate_envelope(envelope)
    # The unique subject constraint plus get_or_create handles first-writer races.
    # A PostgreSQL SELECT FOR UPDATE serializes subsequent reviews for this subject.
    ledger, _ = PreparationLedger.objects.get_or_create(subject_key=envelope['accession_id'])
    ledger = PreparationLedger.objects.select_for_update().get(pk=ledger.pk)
    revision, sequence = envelope['revision'], envelope['review_sequence']
    historical = PreparationOutbox.objects.filter(ledger=ledger, revision=revision, review_sequence=sequence).first()
    if historical:
        return historical
    if sequence != ledger.review_sequence + 1 or (envelope['expected_previous_revision'] or '') != ledger.revision:
        raise PreparationConflict('Stale or out-of-order metadata review')
    project_guid, token = apply_metadata(ledger)
    if not project_guid or not isinstance(token, str) or not re.fullmatch('[a-f0-9]{64}', token):
        raise PreparationConflict('Metadata adapter did not return a valid state token')
    ledger.review_sequence = sequence
    ledger.revision = revision
    ledger.request_id = envelope['request_id']
    ledger.project_guid = project_guid
    ledger.metadata_state_token = token
    ledger.save()
    payload = {
        'request_id': envelope['request_id'], 'revision': revision,
        'review_sequence': sequence, 'accession_id': ledger.subject_key,
        'project_guid': project_guid, 'status': 'prepared', 'prepared': True,
        'searchable': False, 'manual_load_required': True,
    }
    return PreparationOutbox.objects.create(ledger=ledger, revision=revision, review_sequence=sequence, payload=payload)


def deliver_acknowledgment(outbox_id, publish):
    """publish must create an immutable ack, or verify identical bytes on retry.

    A crash after publish but before delivered_at is safe: retry republishes the same
    immutable acknowledgment. This never rolls back clinical state after commit.
    """
    event = PreparationOutbox.objects.get(pk=outbox_id)
    if event.delivered_at is not None:
        return event.payload
    publish(event.payload)
    PreparationOutbox.objects.filter(pk=outbox_id, delivered_at__isnull=True).update(delivered_at=timezone.now())
    return event.payload
