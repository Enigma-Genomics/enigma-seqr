# Offline preparation candidates

family_prep_plan.py validates schema v2 per-person pedigree/HPO/onset/inheritance
without DB access. Source observations must be independently established from pinned
objects and explicit dataset provenance, not copied from the manifest. It emits a
plan only. It does not perform family merge, create individuals or load data.

normalize_singleton_vcf.py requires pysam==0.23.3. Run as a module with source,
NEW output directory, exact source sample ID, and target sample ID. Uses streaming
lossless text copying into BGZF, fresh tabix indexing, sequential and indexed record
count checks, and SHA256 receipts. No original indexes are trusted. The tool accepts
only an explicitly allowlisted GRCh38 reference header with contig checks and final SNV/indel VCFs. Missing, phased,
multiallelic and haploid genotypes stay unchanged; unsupported polyploid, symbolic,
gVCF and unsorted inputs fail without a validation receipt. Failed directories are
not safe outputs; only a completed validation.json authorizes downstream review.
No automatic upload or clinical writes occur.

Legacy Hail SNV_INDEL entries still retain GT/AD/GQ only; DP/PL/PS remain absent from
search exports. Keeping them in staged VCFs does not extend ClickHouse schema.
Merge of final VCFs is NOT joint genotyping. An absent locus must remain a missing
call and cannot establish homozygous reference or reliable de novo status.
Preserve manual Load Data and serial annotation. Normalized, prepared, loaded and
searchable are separate states. Multi-sample Exomiser, kinship/sex QC, resource sizing
and deployed consumer compatibility are release gates, not validated capabilities.


Current reviewed candidate limits: source <=2,000,000,000 bytes, expanded VCF
<=512 MiB, each line <=4 MiB; worker normalization subprocess deadline 1800 seconds.
Worker requires source-size plus 1.1GB temporary disk headroom. These bounds can
reject WGS inputs; no production-sized compatibility or cost claim is established.
Reference headers accept only GRCh38/hg38 and explicit GRCh38.fa/GRCh38.fasta,
hg38.fa/hg38.fasta/Homo_sapiens_assembly38.fasta basenames, plus a declared chr1/1
length of 248956422. This detects obvious mismatches, not independent sequence
provenance; trusted pinned-source build evidence remains required.
