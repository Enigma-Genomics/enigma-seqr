# Disabled family/joint-callset adapter

This is a separate schema, namespace and adapter, not a singleton envelope with
three samples appended. It only represents one already jointly genotyped input.
No clinical activation, migrations, source reads, upload, load trigger, dictionary
refresh or report release is included. `PREPARATION_FAMILY_ENABLED` defaults false
when absent. Test-only overrides do not enable the worker or deployment.

## Exact identity rule (Python and JavaScript)

The server-authoritative joint report ID must match
`[A-Za-z0-9][A-Za-z0-9_.-]{0,127}`. Hash its exact UTF-8 bytes with SHA-256; retain
all 64 lowercase hexadecimal digits. No trimming, case folding, Unicode
normalization or JSON wrapping is permitted.

- Family ID and owned project name: `JR_` + digest (67 characters).
- Ledger subject key: `family:` + digest (71 characters).
- Raw joint report ID is retained separately; it is not an accession ID.

The actual Django Family family_id limit is **100**, not 128; this rule fits both.
The colon cannot occur in accepted singleton accession IDs, keeping ledger keys
separate. Source IDs, target sample IDs and per-person accession IDs use the same
ASCII128 identifier constraint for this family contract. This is intentionally
narrower than arbitrary VCF sample names; unsupported IDs must be mapped/reviewed
explicitly, never silently truncated.

## Envelope and integration API

`preparation_family_contract.validate_family_envelope(envelope)` returns the
validated namespaced subject key. `family_request(envelope)` returns immutable
request content. Use the existing recursive sorted-key `canonical_digest` for
its request_id and the complete envelope excluding revision for revision.
`family_id` and `family_subject_key` implement the shared identity rule.

Exact complete synthetic fixture: `preparation_family_test_utils.synthetic_family_envelope()`.
The schema is version 3, kind `family_joint_callset`; schema 2 singleton behavior
must remain unchanged. Fields include:

- joint_report_id, family_id, subject_key, intake_request_id, request_id, revision;
  intake_request_id identifies the pre-normalization queue, while request_id binds
  the resulting normalized output generations and the intake identity.
- review_sequence and expected_previous_revision with existing ordered CAS rules;
- genome_version `38`, sample_type WES/WGS, dataset_type SNV_INDEL;
- source bucket/path/generation and source_sha256;
- normalized_outputs.vcf and .index, each with exact bucket/path/generation/sha256;
- vcf_sha256/index_sha256 matching those output pins, plus records;
- joint_genotyping_attestation: confirmed=true, authority user_confirmed or
  trusted_workflow, confirmed_by, reference;
- metadata_review: confirmed=true, reviewer, reference. This is a separate explicit
  review of all three people's sex/affected/HPO values, not inferred from calling;
- exactly three individuals with explicit role, sample_id, accession_id,
  source_sample_id, mother_id, father_id, sex, affected and metadata;
- metadata contains exactly hpo_present and hpo_absent arrays. Other clinical
  fields are deliberately unsupported until their mapping is separately reviewed;
- status awaiting_consumer, prepared=false, searchable=false,
  manual_load_required=true, clinical_ready=false,
  calling_history_independently_verified=false.

Roles are proband/mother/father, independent of input order and sex. Only the
proband has parent links; they point to the explicit mother/father sample IDs.
Request fingerprinting orders members proband, mother, father and excludes
reviewable sex/affected/HPO values. It freezes source, normalized output pins,
member/accession/source-sample/parent identities, roles, assay/build, and
joint-genotyping attestation. Replacing any of these needs a separately designed
re-preparation path, not a metadata revision. Hashing/attestation validation does
not independently authenticate clinical review or calling history.

The generic ledger owner must dispatch family validation before selecting the
namespaced subject. The common transaction should retain locking, CAS, exact
retry/outbox replay semantics, and immutable request checks. Family acknowledgments
must carry kind, joint_report_id, family_id and subject_key rather than pretending
the family subject is a patient accession. Backend bridging must use this exact
schema and authoritative report ID; no family cloud producer/consumer is enabled.

## Atomic adapter and ownership

`apply_family_metadata(envelope, ledger, user)` is a transaction-bound callback
for that locked ledger and returns (project_guid, metadata_state_token). It
creates an owned project, one family, three individuals and explicit parent links
atomically. It uses standard application model creation/audit helpers; no existing
project permission or source-storage permission is changed. Standard new-project
owner groups created by Project.save are part of project creation, and need review
before any future activation. No current clinical records are touched by this work.

Existing projects cannot be adopted. Subsequent revisions require matching owned
project GUID/build, exact one-family/three-person membership, unchanged immutable
request ID, and unchanged snapshot token covering every concrete project/family/
individual field (including parent links). Active or inactive datasets block
updates. HPO clearing is explicit. Errors roll back the entire callback and,
through the generic outer transaction, ledger/outbox updates. The gate defaults
off; this does not replace scientific or clinical approval.

Database row locks do not by themselves prove exclusion against every external
loader/manual-edit path. Concurrency with live dataset loading and all clinical
writers remains an integration gate. No loaded-family refresh is supported.

## Verification and remaining gates

Run pure synthetic contract tests in an isolated process:

```sh
PYTHONPATH=. python seqr/tests/standalone/preparation_family_contract_test.py
```

Run actual model tests only through the normal Django **disposable test database**
after generic family dispatch integration:

```sh
python manage.py test seqr.utils.preparation_family_tests
```

Contract tests currently pass (5 cases); syntax and Ruff 0.11.3 checks pass.
Actual-model tests are authored for 10 cases (creation/parent links, idempotence,
HPO clear, new/revision rollback, adoption denial, external edits, loaded data,
extra family, immutable assay, disabled gate), but local minimal environments do
not contain full application dependencies and they have not been run here.
They must pass the application CI/disposable PostgreSQL environment before this
adapter can be considered verified. Generic ledger tests are owned separately.

Further gates: independent audit; exact frontend/backend schema and golden vectors;
authoritative report ownership; pinned source/output verification by the future
consumer; real provenance, sex/kinship and metadata review; multi-sample Exomiser;
concurrent writer/loader coordination; permissions review; production resource
supervision; and explicit authorization for any live clinical activation/load.

The singleton validator reserves the exact `JR_[a-f0-9]{64}` project namespace.
Manual external project creation is not serialized by this ledger and remains a
runtime activation gate. The adapter's atomic decorator does not itself prove
that the caller holds the ledger lock; callers must use `apply_preparation`.
