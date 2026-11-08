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

Default wall timeout: 1800 seconds. Limits: 2 GB compressed/plain source bytes,
1 million records per source, 512 MiB combined expanded text, 1 MiB per line,
32 alleles per source site. Validation stores records in memory, not constant-memory
streaming. These are work bounds, not a hard RSS guarantee; use an appropriately
memory/disk-limited container before service integration. Current expanded limit
may reject the proposed 3 x 150 MB compressed inputs. Do not raise it for production
without representative synthetic scale testing and resource sizing.

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
