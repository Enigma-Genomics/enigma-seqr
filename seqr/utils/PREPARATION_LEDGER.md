# Versioned preparation ledger candidate

Not deployed. Migration 0091 has not been executed against an operational database.
Worker activation remains blocked on end-to-end authority and runtime validation.

An immutable envelope has schema_version=2, kind=individual, a content-hash revision,
review_sequence (first 1, thereafter +1), and expected_previous_revision (null for
first review, otherwise the last applied hash). The revision is SHA256 of UTF-8
canonical JSON with sorted object keys, compact separators, excluding revision.
The producer must derive order from a trusted server-side review sequencer, not a
client clock or lexicographic hashes. No source currently establishes that sequencer.

PreparationLedger uniquely keys the accession, locks it with SELECT FOR UPDATE,
compares the expected predecessor, and commits the metadata adapter, new revision,
and PreparationOutbox acknowledgment in ONE Postgres transaction. Duplicate applied
revisions replay their original outbox event without reapplying clinical metadata.
Historical acks include their sequence; consumers must never treat an older ack as
the newest review's success. A future/gapped or conflicting revision changes nothing.

The restricted individual adapter locks the project and individual, detects edits
outside the ledger via a full metadata/last-modified fingerprint, and rejects initial
adoption of existing clinical rows or updates to already-loaded individuals. A
separately reviewed baseline-adoption and loaded-data refresh procedure is required
before those cases can be enabled. Candidate-gene free text is not silently stored
in a structured field: nonempty free-text candidate genes are rejected pending a
validated ID mapping. Missing optional metadata preserves the prior value; explicit
empty HPO arrays clear HPO lists. Unknown sex/affected remain unknown.

Dictionary refresh is deliberately outside the clinical transaction. After commit,
the ack dispatcher refreshes dictionaries and creates an immutable cloud ack. It
marks the outbox delivered only after success. Crashes after external publication
are replayable because publication must create-if-absent or verify identical bytes.
Prepared means metadata transaction plus dictionary refresh acknowledged; it NEVER
means loaded/searchable. Manual Load Data and serial annotation remain unchanged.

Standalone tests use synthetic tables only, not migrations or clinical rows:
PYTHONPATH=. python seqr/tests/standalone/preparation_ledger_test.py
SQLite covers rollback/idempotency/replay; it does NOT prove row-lock behavior.
For ephemeral CI Postgres only, PREP_TEST_POSTGRES_PORT connects exclusively to
127.0.0.1, database prep_ledger_synthetic, user/password synthetic. This creates only
synthetic test tables through schema_editor. Never point tests at operational DBs.
No arbitrary DSN/environment credential loading is supported.

Ordinary metadata revisions pin the original normalization request ID, which binds
accession, assay, reference and immutable source identity/header sample. A changed
file or generation is NOT an ordinary retry: it is rejected even before loading.
An explicit re-preparation workflow with prior-source/revision checks and a verified
no-loaded-data condition remains separate and disabled. No implicit source switch
is permitted merely because an accession matches.

Validation checkpoint: the initial migration and three actual Individual adapter
tests passed in disposable application CI databases; PostgreSQL standalone CAS/race
tests passed separately. No production migration or activation has occurred.
