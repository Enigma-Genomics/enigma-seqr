# Experimental offline final-VCF trio union

`merge_trio_vcfs.py` is an opt-in preparation primitive, not a production loading
path or joint-genotyping workflow. Requires the vetted `pysam==0.23.3` installation
with bundled bcftools. No network, patient records, uploads, or SEQR mutations.

Run with the helper's Python environment:

```
python -m loading_pipeline.ops.merge_trio_vcfs local-manifest.json NEW_OUTPUT_DIR
```

Local manifest (separate from family-plan schema v2; the caller must validate and
translate the approved family manifest):

```
{
  "genome_version": "GRCh38",
  "members": [
    {"role": "proband", "path": "/local/child.vcf.gz", "vcf_sample_id": "original_child", "sample_id": "child", "genome_version": "GRCh38"},
    {"role": "mother", "path": "/local/mother.vcf.gz", "vcf_sample_id": "original_mother", "sample_id": "mother", "genome_version": "GRCh38"},
    {"role": "father", "path": "/local/father.vcf.gz", "vcf_sample_id": "original_father", "sample_id": "father", "genome_version": "GRCh38"}
  ]
}
```

Each member may supply `expected_sha256`, checked before processing. Include
bucket/path/generation provenance in each member for the receipt. This offline
helper cannot verify cloud generations: the download/orchestration layer must pin
and verify those itself. Local bytes are SHA256 hashed before and after processing.

Input contract is intentionally narrow and fail-closed:

- Three final singleton VCFs; explicit roles and distinct target IDs/paths.
- Exact `##reference=GRCh38` or `GRCh37` matching manifest, identical ordered
  contig names and lengths. A FASTA URL or filename is NOT inferred as a build.
  Obtain verified build metadata upstream; do not rewrite it based on a guess.
- Sorted, at most one record per source locus. No structural/symbolic alleles,
  gVCF `<NON_REF>` / `<*>`, END blocks, or ambiguous duplicate loci. gVCFs require
  a separately validated joint-genotyping workflow.
- GT is identified by FORMAT name; reordered GT is explicitly moved first for
  standards-compliant bcftools parsing. Haploid/diploid, partial-missing and phased
  GT are supported. AD and PL vectors are verified after allele remapping.
- **Allele-dependent INFO Number=A/R/G is currently rejected.** Shared sites must
  agree exactly on REF/ID/QUAL/FILTER/scalar INFO. Common GATK outputs contain AC/AF
  or conflicting site annotations and will fail here. Never strip those fields
  to force a pass: extending reconciliation requires separate tests and review.
- Source FORMAT values and every source locus are checked against merged output;
  no variant or unsupported field is silently dropped. Missing sites remain missing
  and are never converted to homozygous reference. This is a union, not re-calling.

Output is a new exclusively created directory with `callset.vcf.gz`, fresh
`callset.vcf.gz.tbi`, and `validation.json` written last. Existing destinations are
not overwritten; exceptions/timeouts clean the newly owned output and scratch
files. Do not consume output before successful return and receipt verification.
Input indexes are ignored; indexes for temporary copies and merged output are
regenerated and output indexed-read counts checked. OS termination or disk failure
may leave incomplete files; absence of validation.json means unusable output.

Default wall timeout: 1800 seconds. Default bounds: 2 GB total compressed/plain
source bytes, 20 million records per source, 8 GiB combined expanded text, 8 GiB
scratch files, 1 MiB per line, 16 MiB header per source, 32 alleles per source site.
`--max-records`, `--max-expanded-bytes`, and `--max-staging-bytes` configure the
bounds; expanded/scratch configuration cannot exceed 64 GiB. They are explicit
resource allowances, not a promise that the machine has that capacity.

Validation now merge-joins three sorted singleton readers against the merged
reader, retaining one current record per source rather than all records/loci.
Headers and rows have explicit size bounds. The supervisor enforces wall timeout,
polls aggregate scratch size and free space every 50 ms, kills and reaps the
worker on breach, and uses the OS per-file size limit as an additional guard.
Aggregate usage can briefly overshoot by writes within one polling interval;
use a filesystem/container quota for a strict total-storage ceiling. The free
space floor is 256 MiB. Input source files are not charged to the scratch budget.
The receipt records configured limits and the sampled scratch peak.

A 450 MB compressed input set fits the compressed-byte default, but expansion and
number of rows must still fit the configured bounds. Raising limits is not a
scientific compatibility approval. Cloud Run writable temporary storage may count
against container memory: deployment must provide and size appropriate storage
and memory; these local tests do not validate a Cloud Run configuration.

Receipt explicitly states `joint_genotyped=false`, `prepared=false`,
`searchable=false`, manual loading required. No kinship, inferred sex, pedigree
biological validation or clinical validation is performed. Exomiser trio remains
unsupported until separately tested end-to-end.

## Checks

```
python -m pytest -q -c /dev/null loading_pipeline/ops/merge_trio_vcfs_test.py
ruff check loading_pipeline/ops/merge_trio_vcfs.py loading_pipeline/ops/merge_trio_vcfs_test.py
```

Synthetic-only local smoke measurement on Bs-M2: 30,000 input rows (three identical
10,000-row, 1,497,789-byte total uncompressed sources) produced 10,000 union rows in
0.686 seconds; child max RSS 58,114,048 bytes. This is not a benchmark of GATK
annotations, 450 MB compressed input, cloud network/upload, SEQR loading, cost or
clinical workloads. No time/cost extrapolation is justified from this small test.


## Scientific annotation policy: decision remains required

The tool uses the pinned bundled `bcftools merge -m all -i -`, without
`--missing-to-ref`, `--gvcf`, a field-removal operation, or default DP/DP4 summation.
The bundled merge help explicitly documents `-i -` as disabling default INFO
aggregation, and missing-vector default `.`. Disabling INFO aggregation does NOT
make bcftools a statistical joint genotyper. Output is checked against every
source record and fails rather than accepting changed site values.

Representative synthetic GATK-like tests now name the actual gate: INFO AC/AF
(Number=A) are rejected with their names; conflicting scalar QD produces an
INFO/QD policy error. The following common fields require a scientific decision:

- AC/AF/AN: recounting from observed GT would describe only the available calls,
  not re-establish uncalled-reference evidence, original cohort counts or joint
  calls. Missing sites must never increase the denominator as presumed 0/0.
- MLEAC/MLEAF and allele-specific AS_* vectors: final GT alone cannot reconstruct
  the callers' likelihood estimates or underlying allele-specific evidence.
- QUAL, QD, FS, SOR, MQ, MQRankSum, ReadPosRankSum, BaseQRankSum, and raw MQ metrics:
  selecting one value, taking a maximum, averaging or summing is not a generally
  valid combined statistic. Required sufficient statistics may not be present.
- VQSLOD/AS_VQSLOD, culprit and FILTER labels: model/filter provenance is specific
  to the original run. Union, intersection, selecting PASS or setting FILTER to
  missing can change which variants downstream considers usable.
- INFO DP: summation is an explicit policy choice, not permission implied by
  `merge`. Sample FORMAT DP/AD/GQ/PL stay per sample and are verified separately.

**Current decision:** retain the strict gate. Do not strip INFO fields or unset
QUAL/FILTER to force a pass. Biswa/scientific owners must choose validated
joint-genotyped input or approve a specific final-VCF annotation interpretation
with downstream SEQR/Exomiser validation before clinical use.

A possible future proposal is to retain exact source records in a hashed,
sample-keyed provenance sidecar and avoid claiming joint site statistics. This
is NOT implemented or approved. It requires explicit review of downstream missing
FILTER/QUAL/INFO handling, access/provenance retention, and clinical presentation.
The current work does not weaken the gate or activate that alternative.

## Reproducible bounded scale smoke test

```
python -m loading_pipeline.ops.benchmark_trio_merge --rows 10000
python -m loading_pipeline.ops.benchmark_trio_merge --rows 200000
python -m loading_pipeline.ops.benchmark_trio_merge --rows 350000
```

The generator writes only synthetic BGZF files, bounded to at most one million
rows per source, with a deterministic high-entropy scalar annotation. All loci
are shared, annotations agree, and genotypes have AD/PL; this stresses streaming
and compression, not clinical annotation compatibility. Generation time is
excluded from merge wall time. Temporary inputs and outputs are deleted afterward.

On Bs-M2, the streaming implementation measured:
- 10,000 rows/source: 9,151,938 compressed / 17,127,921 expanded input bytes;
  1.720 seconds merge wall time, 36,388,864 bytes child max RSS.
- 200,000 rows/source: 183,062,940 compressed / 343,467,924 expanded input bytes;
  29.773 seconds merge wall time, 33,374,208 bytes child max RSS.

- 350,000 rows/source: 320,351,586 compressed / 601,317,924 expanded input bytes;
  54.550 seconds merge wall time, 33,538,048 bytes child max RSS, sampled scratch
  peak 427,544,410 bytes. This exceeds the former 512 MiB expanded-text cap.

These local measurements show record-count-independent validation memory in
these fixtures. They do not predict GATK policy compatibility, disparate alleles,
cloud vCPU throughput, network/download time, SEQR ingestion or monetary cost.
