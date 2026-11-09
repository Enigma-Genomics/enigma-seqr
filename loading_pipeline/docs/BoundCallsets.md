# Bound callset inputs

A loading request may include `callset_identity` with `vcf` and `index` objects.
Each object contains `path`, `generation`, `size`, and `crc32c`; an optional
publisher-attested `sha256` is accepted. Generation and size are positive decimal
strings, CRC32C is canonical base64, and SHA256 is lowercase hexadecimal. The
index URI must equal the VCF URI plus `.tbi`. Legacy single-file VCFs may use
`index: null`; managed aliases must have an indexed pair. Bound inputs must be exact GCS
objects, not wildcard shards or MatrixTables.

The application resolves a managed fixed-name alias to its unique prepared VCF
and index and sends those identities. The pipeline does not import the mutable
alias. GCS metadata is checked with `ifGenerationMatch` before enqueue validation,
before and after task completion checks, and before and after task execution.
For Hail imports, the execution check surrounds the materializing write, not just
the lazy import expression. Dynamic Luigi tasks are checked again on resumption.
Missing objects and failed metadata preconditions fail closed. SHA256 is an
attestation from the publisher, not independently recomputed by the pipeline;
generation, CRC32C, and size are compared against GCS on each check.

All callset-level cache paths include the complete canonical VCF/index identity.
Replacing either object, including a replacement of equal size, cannot reuse the
old identity's cache. Queue retries retain the identity unchanged.

Set `CALLSET_BINDING_ENABLED=true` (or `1`) to reject unbound managed aliases.
An unbound single GCS VCF, or each member of an unbound GCS wildcard, is checked
for the `seqr_binding` metadata key. Even an empty or malformed value requires
binding. Unrelated legacy VCFs, wildcard shards and MatrixTables remain supported.
Supplied identities are always verified regardless of the flag.

Before any bound task can reuse or write outputs, an immutable run manifest binds
the run ID to the canonical input identity, genome, dataset type, sample type,
projects, skipped validations, and sex/relatedness and TDR flags. Retry attempt
numbers are excluded. GCS uses a create-only generation precondition; local runs
publish a fully written manifest with an atomic non-replacing link. Existing run
outputs without a manifest cannot be adopted by a bound request. A mismatching
manifest rejects run-ID reuse. An uncertain manifest upload is reconciled by
reading its binding before any work continues.

Roll out compatible app and pipeline code before enabling the flag. No new
credentials, IAM grants, or storage configuration are required.
