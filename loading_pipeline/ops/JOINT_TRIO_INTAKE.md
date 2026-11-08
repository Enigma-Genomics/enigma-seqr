# Offline intake of an already jointly genotyped trio

`normalize_joint_trio_vcf.py` is separate from the experimental union of three
final singleton VCFs. It accepts exactly one already jointly genotyped three-sample
VCF. It does not merge, joint genotype, recalculate QUAL/FILTER/INFO, impute absent
calls, sort records, or change any data row. Only the three sample labels in the
column header change; their source order is preserved.

This is **offline, local validation only**. No database, storage bucket, registry,
ledger, report, load trigger, or clinical state is contacted or changed. The
receipt always says `prepared=false`, `searchable=false`, `clinical_ready=false`,
`manual_load_required=true`, and `family_namespace_integrated=false`.

## Inputs and trust boundary

Use Python 3.11 and pysam 0.23.3 (the synthetic verification environment). All
inputs must be local files. A caller must supply an independently obtained,
verified GRCh38 reference dictionary and explicit joint-genotyping attestation.
A header containing GATK or any other tool name is **not** provenance proof.
The helper cannot independently verify either human/workflow attestation or
external evidence references. Its receipt retains those assertions and explicitly
sets `calling_history_independently_verified=false`, `kinship_verified=false`,
and `sex_verified=false`.

Example manifest (synthetic identities only; replace the hash with the SHA-256
of the exact input bytes):

```json
{
  "schema_version": 1,
  "kind": "attested_joint_trio",
  "family_id": "SYNTH-FAMILY",
  "joint_report_id": "SYNTH-REPORT",
  "genome_version": "38",
  "sample_type": "WGS",
  "dataset_type": "SNV_INDEL",
  "source_sha256": "<64 lowercase hexadecimal characters>",
  "joint_genotyping_attestation": {
    "confirmed": true,
    "authority": "user_confirmed",
    "confirmed_by": "synthetic-operator",
    "reference": "synthetic-run-record"
  },
  "individuals": [
    {"role":"proband", "source_sample_id":"CHILD", "sample_id":"SYNTH-P", "accession_id":"SYNTH-AP", "father_id":"SYNTH-F", "mother_id":"SYNTH-M", "sex":"U", "affected":"U", "metadata":{"hpo_present":[], "hpo_absent":[]}},
    {"role":"mother", "source_sample_id":"MOM", "sample_id":"SYNTH-M", "accession_id":"SYNTH-AM", "father_id":null, "mother_id":null, "sex":"U", "affected":"U", "metadata":{"hpo_present":[], "hpo_absent":[]}},
    {"role":"father", "source_sample_id":"DAD", "sample_id":"SYNTH-F", "accession_id":"SYNTH-AF", "father_id":null, "mother_id":null, "sex":"U", "affected":"U", "metadata":{"hpo_present":[], "hpo_absent":[]}}
  ]
}
```

`authority` is `user_confirmed` or `trusted_workflow`. Identity fields are explicit,
unique within each category, case-sensitive ASCII identifiers. Family and report
IDs are separately pinned as CLI arguments and must match the manifest and differ
from each other. This is consistency validation, not allocation/reservation of
an identity in any downstream system. Roles never imply sex or sample order.
Sex accepts M/F/U and affected accepts A/N/U. Metadata is retained per person;
HPO lists are validated for syntax, duplicates, and contradictory entries.
Other metadata is retained as an opaque assertion, not clinically validated.

Example reference evidence for the one-contig synthetic fixture:

```json
{"genome_version":"38", "verified":true, "verification_reference":"synthetic-reference-fixture", "contigs":{"1":248956422}}
```

For a real calling reference, provide its complete sequence dictionary, not this
fixture and not a dictionary copied unchecked from the input. Every header contig
must match name and length exactly and the complete dictionary must match. The
helper accepts explicit GRCh38/hg38 reference aliases and requires the canonical
chromosome 1 length. It does not verify sequence content against a FASTA; verified
reference derivation remains the caller's responsibility. Contig order follows
the source header; records must already be sorted by that order and position.
No chromosome alias conversion is performed.

## CLI and outputs

```sh
python loading_pipeline/ops/normalize_joint_trio_vcf.py \
  source.vcf.gz NEW_OUTPUT_DIRECTORY manifest.json reference.json \
  SYNTH-FAMILY SYNTH-REPORT \
  --timeout 1800 --max-source-bytes 536870912 \
  --max-expanded-bytes 536870912 --max-line-bytes 4194304 \
  --max-header-bytes 8388608 --max-output-bytes 536870912 \
  --max-scratch-bytes 1073741824
```

Plain VCF, gzip and BGZF inputs are accepted. Source indexes are ignored. Output
is `callset.vcf.gz`, a fresh `.tbi`, and `validation.json` written last. Sequential
pysam parsing and index-based enumeration must agree with input record count;
index row content is SHA-256 checked against the original row stream. A test
also checks a bounded indexed region. Input SHA-256 is checked before and after
processing; output VCF/index hashes are recorded. Consumers must not inspect an
in-progress directory as if it were a completed artifact. Existing output
directories are refused; caught failures remove the directory created by this
invocation. Abrupt process/host termination can leave an incomplete directory;
no receipt means no completed result. No restart/resume is supported.

Missing, partially missing, phased and haploid/diploid GTs are preserved, including
nonvariant rows with ALT `.`. Polyploid/malformed/out-of-range GTs, gVCF END blocks,
symbolic alleles, breakends, spanning deletions, missing/extra/duplicate sample IDs,
build/contig mismatch, unsorted positions and unsupported headers fail closed.
No records are silently dropped. This is not full VCF semantic or clinical QC.

## Bounds and verification

The input is streamed with bounded binary line reads. JSON inputs are capped at
1 MiB each. Positive CLI limits cover source/compressed bytes, expanded bytes,
line bytes, header bytes, each output file, aggregate scratch, and walltime.
The helper conservatively reserves four times expanded bytes plus 1 MiB against
both output and scratch limits, so default effective input capacity is smaller
than the expanded limit. This reservation is deliberately conservative, not a
production capacity claim. Indexing has a subprocess timeout. Hashing, parsing,
writing and index enumeration check elapsed time between bounded operations;
this is not an OS-level hard limit on blocked filesystem/native calls. Production
use would additionally need a job supervisor with hard time/disk/memory limits,
validated workload sizing and operational cleanup. No cloud benchmark was run.

Synthetic verification:

```sh
cd loading_pipeline/ops
python -m unittest normalize_joint_trio_vcf_test -v
# from repository root, using Ruff 0.11.3:
ruff check loading_pipeline/ops/normalize_joint_trio_vcf*.py
ruff format --check loading_pipeline/ops/normalize_joint_trio_vcf*.py
```

Tests include 20,000 synthetic records, lossless annotations/GT preservation,
source gzip, fresh index despite a stale source index, region retrieval, contract
and malformed-record rejection, resource limits, index-timeout cleanup, output
non-overwrite and CLI invocation. No patient data or cloud access is used.

## Remaining gates

Independent code/security/scientific review and owner integration are required.
This helper does not solve the separate final-singleton union annotation policy.
Actual calling history and reference evidence must be confirmed separately.
Kinship, sex QC, clinical metadata review, family/report/sample namespace design,
consumer contracts, multisample Exomiser behavior, representative production
capacity, and explicit load/release authorization remain unvalidated. Do not wire
this receipt into singleton ledgers or consumers. No live/deployment activation
is included or implied.
