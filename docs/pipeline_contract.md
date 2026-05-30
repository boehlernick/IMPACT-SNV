# IMPACT-SNV Pipeline Contract

## Purpose

This document defines the current IMPACT-SNV pipeline contract before the FAVOR-CLI refactor.

It describes the existing workflow, inputs, outputs, file naming conventions, and behavior that should be preserved unless a future refactor explicitly changes it.

This file is intended to support the staged FAVOR-CLI refactor by giving developers and GitHub Copilot coding agents a stable description of the current pipeline.

This document should be read alongside:

- `.github/copilot-instructions.md`
- `docs/favorcli_refactor_plan.md`
- `docs/annotation_compatibility_contract.md`
- `docs/output_gds_contract.md`
- `docs/chromosome_handling_contract.md`

## High-Level Workflow

The current IMPACT-SNV workflow is organized as four sequential steps:

```text
Input VCF / VCF.GZ files
→ Step 1: VCF merge, normalization, and chromosome split
→ Step 2: VCF-to-GDS conversion
→ Step 3: backend-selected annotation (legacy FAVORannotator compatibility or FAVOR-CLI path)
→ Step 4: IMPACT variant prioritization
→ Final per-sample *_SNV_IMPACT.gds files
```

Current step folders:

```text
step1_vcf_merge/
step2_vcf2gds/
Step3_favorannotator-rap/
step3_favorcli_annotation/
impact_snv/favor/
step4_impact_prioritization/
```

The current workflow is primarily implemented as DNAnexus/RAP applets, but the refactor should progressively separate core logic from platform wrappers where practical.

## Current Pipeline Contract Summary

### Step 1: VCF Merge, Normalization, and Chromosome Split

Folder:

```text
step1_vcf_merge/
```

Primary script:

```text
step1_vcf_merge/src/IMPACT_SNV_indel.py
```

Current role:

- Accept one or more VCF or VCF.GZ input files.
- Optionally accept VCF index files.
- Accept a reference genome FASTA for normalization.
- Download inputs in a DNAnexus/RAP execution environment.
- Ensure inputs are bgzipped.
- Index VCF files when indexes are not provided.
- Check VCF headers with `bcftools view`.
- Remove FORMAT fields other than `GT`.
- Merge input VCF files.
- Normalize and split multiallelic variants using `bcftools norm`.
- Split the normalized merged VCF into chromosome-specific VCF.GZ files.

Current important behavior:

- Input genotype data are reduced to `GT` FORMAT only.
- The merged output is bgzipped and indexed.
- Chromosome-specific VCF files are emitted using the `merged_chr*.vcf.gz` naming convention.

Current chromosome split behavior:

- The split step attempts chromosome-specific outputs for `chr1` through `chr22`, `chrX`, and `chrY`.
- The region naming used for `bcftools view --regions` matches the chromosome naming style already present in the merged VCF.
- Chromosomes with no variants are skipped with a log message instead of failing the step.

Current remaining limitations:

- Mitochondrial chromosomes are not part of the current pipeline contract.
- Numeric sex chromosome aliases such as `23` and `24` are not silently converted.

Current Step 1 input contract:

```text
input_vcfs: one or more .vcf or .vcf.gz files
index_files: optional .tbi or .csi files
reference_genome: .fa or .fasta file
```

Current Step 1 output contract:

```text
merged_vcf: merged_output.vcf.gz
split_vcfs: array of merged_chr*.vcf.gz files
```

Expected current output naming pattern:

```text
merged_chr1.vcf.gz
merged_chr2.vcf.gz
...
merged_chr22.vcf.gz
merged_chrX.vcf.gz
merged_chrY.vcf.gz
```

## Step 2: VCF-to-GDS Conversion

Folder:

```text
step2_vcf2gds/
```

Primary scripts:

```text
step2_vcf2gds/src/vcf2gds.sh
step2_vcf2gds/resources/home/dnanexus/vcf2gds.R
```

Current role:

- Accept one or more chromosome-specific VCF or VCF.GZ files.
- Convert each VCF file to SeqArray GDS format using `SeqArray::seqVCF2GDS`.
- Add `annotation/info/QC_label` with value `PASS` for every variant.
- Upload generated GDS files as DNAnexus applet outputs.

The current input/output naming path is compatible with X/Y chromosome files emitted by Step 1.

Current important behavior:

- The VCF-to-GDS conversion is performed by R using SeqArray.
- The GDS file becomes the core internal exchange format for downstream annotation and prioritization.
- The `QC_label` node is added under `annotation/info` for all variants.

Current Step 2 input contract:

```text
vcf_files: one or more .vcf or .vcf.gz files
```

Current Step 2 output contract:

```text
gds_files: array of .gds files
```

Expected current output naming pattern:

```text
merged_chr1.gds
merged_chr2.gds
...
merged_chr22.gds
```

Target refactor note:

Step 2 should remain mostly stable during the initial FAVOR-CLI refactor. If X/Y VCF files are produced by Step 1, Step 2 should process them through the same conversion path where possible.

## Step 3: Annotation Backends

Current implementation includes a backend selector in:

```text
impact_snv/favor/annotate.py
```

exposed by:

```text
impact-snv favor-annotate --backend ...
```

Current backend values:

```text
favor-cli
legacy-favorannotator
favor-cli-skeleton
```

The default remains `favor-cli`.

### Step 3A: Legacy FAVORannotator Compatibility Backend

Folder:

```text
Step3_favorannotator-rap/
```

Primary scripts:

```text
Step3_favorannotator-rap/src/code.sh
Step3_favorannotator-rap/resources/home/dnanexus/favorannotator.R
```

Current role:

- Accept a chromosome-level GDS file.
- Accept an output prefix.
- Accept a chromosome number.
- Optionally accept a FAVOR CSV database file.
- Run a legacy FAVORannotator-style R workflow.
- Generate variant identity information from the GDS.
- Join GDS variants against FAVOR CSV annotation chunks.
- Add FAVOR-style annotations into the GDS under:

```text
annotation/info/FunctionalAnnotation/
```

Current chromosome/file discovery behavior:

- The filename pattern `^merged_chr.*\.gds$` matches autosomes and sex chromosomes such as `merged_chrX.gds` and `merged_chrY.gds`.
- The existing export path preserves sample-specific files such as `{sample_id}_chrX.gds` and `{sample_id}_chrY.gds` when those chromosomes are present.

Current important behavior:

- Step 3 mutates or rewrites the GDS so that it contains a `FunctionalAnnotation` annotation namespace.
- Step 4 currently depends on this namespace.
- The legacy Step 3 implementation includes assumptions about FAVOR database files, FAVOR CSV schema, xsv-based joins, and chromosome-specific annotation chunks.

Current Step 3 input contract:

```text
outfile: output file prefix
gds_file: chromosome-level .gds file
chromosome: chromosome number
use_compression: YES or NO
favor_csv_database: optional FAVOR CSV database .txt file
```

Current Step 3 output contract:

```text
results: annotated .gds file
```

Expected annotation namespace produced by Step 3:

```text
annotation/info/FunctionalAnnotation/
```

Required Step 4 compatibility fields currently expected under this namespace include:

```text
annotation/info/FunctionalAnnotation/clnsig
annotation/info/FunctionalAnnotation/genecode_comprehensive_info
annotation/info/FunctionalAnnotation/genecode_comprehensive_exonic_category
annotation/info/FunctionalAnnotation/refseq_exonic_category
annotation/info/FunctionalAnnotation/ucsc_exonic_category
annotation/info/FunctionalAnnotation/apc_protein_function_v3
```

Current compatibility note:

The legacy Step 3 should not be deleted or rewritten at this stage. The compatibility backend wraps legacy outputs and stages them into canonical paths for downstream build/finalize steps.
It does not imply that `impact-snv favor-annotate --backend legacy-favorannotator` locally executes the DNAnexus FAVORannotator applet.

### Step 3B: FAVOR-CLI Path

Current FAVOR-CLI-related modules:

```text
impact_snv/favor/ingest.py
impact_snv/favor/annotate.py
step3_favorcli_annotation/
```

`favor-cli-skeleton` currently supports dry-run integration only and does not claim Step 4-ready scientific equivalence.

Native FAVOR `.cohort` or FAVOR-generated genotype outputs remain experimental/discovery-only and are not part of the v1.0.0 release gate. The supported release path is FAVOR annotation-only plus `impact-snv extract-genotypes`, followed by `build-gds`, `finalize-gds`, and `validate-gds`.

Target refactor note:

A production FAVOR-CLI-backed annotation adapter should continue to evolve under:

```text
step3_favorcli_annotation/
```

The FAVOR-CLI path must produce outputs that satisfy the same Step 4 annotation compatibility contract before it can replace legacy production usage.

## Step 4: IMPACT Variant Prioritization

Folder:

```text
step4_impact_prioritization/
```

Primary scripts:

```text
step4_impact_prioritization/src/code.sh
step4_impact_prioritization/resources/home/dnanexus/IMPACT-prioritization.r
step4_impact_prioritization/resources/home/dnanexus/seqarray_append.R
```

Current role:

- Accept one or more annotated chromosome-level GDS files.
- Accept a phenotype-specific gene-disease association file.
- Read required FAVOR-style annotation nodes from each GDS.
- Compute IMPACT prioritization scores.
- Export per-sample chromosome-level GDS files containing prioritized variants.
- Merge per-sample chromosome GDS files into final per-sample files.
- Post-process final files to add IMPACT output annotations.

Current Step 4 input contract:

```text
genelist: GeneList.txt or equivalent tab-delimited file
gds_files: one or more annotated .gds files
```

The gene-disease association file must contain:

```text
symbol
globalScore
```

Current Step 4 required annotation input nodes:

```text
annotation/info/FunctionalAnnotation/clnsig
annotation/info/FunctionalAnnotation/genecode_comprehensive_info
annotation/info/FunctionalAnnotation/genecode_comprehensive_exonic_category
annotation/info/FunctionalAnnotation/refseq_exonic_category
annotation/info/FunctionalAnnotation/ucsc_exonic_category
annotation/info/FunctionalAnnotation/apc_protein_function_v3
```

Current scoring formulas:

```text
Tier 1:
80 + 20 * globalScore

Tier 2:
60 + 40 * globalScore

Tier 3:
20 + 80 * globalScore

Tier 4:
100 * ((0.5 * normalized_APC) + (0.5 * globalScore))
```

These formulas are part of the current scientific contract and must not be changed during the FAVOR-CLI adapter refactor unless explicitly requested.

Current Step 4 output contract:

```text
{sample_id}_SNV_IMPACT.gds
impact_prioritization.log
```

Final output files should contain, at minimum:

```text
annotation/info/impact_score
annotation/info/impact_score_calc
annotation/info/tier
annotation/info/clnsig_flags/*
```

Current output naming pattern:

```text
{sample_id}_SNV_IMPACT.gds
```

Target refactor note:

Step 4 should remain mostly stable until the FAVOR-CLI annotation adapter and GDS injection layer are validated. Minimal future Step 4 changes may include:

- Better preflight validation of required annotation nodes.
- Better error messages for missing or malformed annotations.
- More robust ClinVar token parsing while preserving existing behavior.
- Better X/Y file discovery and merging.

## Platform Contract

The current pipeline is primarily implemented as DNAnexus/RAP applets.

Current platform assumptions include:

- Inputs are often downloaded using `dx download`.
- Outputs are registered using `dx-upload-all-outputs`.
- Each step contains a `dxapp.json` applet specification.
- Runtime environments are specified through DNAnexus applet configuration.

Target refactor direction:

- Preserve DNAnexus/RAP support.
- Avoid making DNAnexus the only execution model for new code.
- Prefer separating reusable core logic from platform-specific wrappers.
- Add local, container, or HPC entrypoints only after the FAVOR-CLI adapter path is stable enough to justify portability work.

## File Naming Contract

Current expected intermediate naming conventions include:

```text
merged_chr1.vcf.gz
merged_chr2.vcf.gz
...
merged_chr22.vcf.gz

merged_chr1.gds
merged_chr2.gds
...
merged_chr22.gds
```

Current final naming convention:

```text
{sample_id}_SNV_IMPACT.gds
```

Planned X/Y-compatible naming conventions:

```text
merged_chrX.vcf.gz
merged_chrY.vcf.gz
merged_chrX.gds
merged_chrY.gds
{sample_id}_chrX.gds
{sample_id}_chrY.gds
```

Changing file naming conventions should be treated as a compatibility-impacting change and should be avoided unless explicitly justified.

## Current Known Limitations

The current pipeline has the following known limitations relevant to the FAVOR-CLI refactor:

1. Step 1 currently splits only autosomes.
2. X and Y chromosome support is incomplete.
3. Legacy Step 3 is tightly coupled to FAVORannotator assumptions.
4. Legacy Step 3 contains static/hard-coded FAVOR database behavior.
5. Step 4 directly depends on legacy-style `annotation/info/FunctionalAnnotation/*` nodes.
6. The repo is strongly DNAnexus/RAP-oriented.
7. There is limited automated validation for schema compatibility.
8. FAVOR-CLI output schema compatibility has not yet been proven.

## Refactor Compatibility Requirements

During the FAVOR-CLI refactor, the following must remain true unless explicitly changed:

- Input VCF expectations remain broadly compatible.
- GDS remains the downstream internal exchange format.
- Step 4 scoring formulas remain unchanged.
- `GeneList.txt` format remains unchanged.
- Final output names remain `{sample_id}_SNV_IMPACT.gds`.
- Final output nodes remain compatible with IMPACT-VIS expectations.
- Legacy FAVORannotator-produced GDS files remain usable by Step 4.

## What Future FAVOR-CLI Work Must Preserve

Any future FAVOR-CLI annotation path must eventually produce GDS files that satisfy the Step 4 input contract.

At minimum, before Step 4 runs, an annotated GDS must contain:

```text
annotation/info/FunctionalAnnotation/clnsig
annotation/info/FunctionalAnnotation/genecode_comprehensive_info
annotation/info/FunctionalAnnotation/genecode_comprehensive_exonic_category
annotation/info/FunctionalAnnotation/refseq_exonic_category
annotation/info/FunctionalAnnotation/ucsc_exonic_category
annotation/info/FunctionalAnnotation/apc_protein_function_v3
```

The values must be:

- Aligned to the GDS variant order.
- The same length as the GDS variant count.
- Readable with `seqGetData()`.
- Typed consistently with Step 4 expectations.
- Validated before prioritization.

## Validation Expectations

Future changes to pipeline behavior should include at least one of:

- Unit tests.
- Fixture-based tests.
- Dry-run validation.
- Preflight schema validation.
- Documented manual validation commands.

For this pipeline contract specifically, validation should eventually confirm:

- Step 1 does not drop expected chromosomes.
- Step 2 preserves sample and genotype information.
- Step 3 or its FAVOR-CLI replacement writes required annotation fields.
- Step 4 produces expected per-sample output files.
- Final GDS files contain expected IMPACT annotations.
- Legacy and adapter-injected GDS files are both supported during the transition period.

## Initial Implementation Guidance for Coding Agents

When using a coding agent, the first implementation task should be documentation-only.

Do not begin by replacing FAVORannotator.

Recommended first coding-agent scope:

```text
Implement Milestone 1 only: documentation and contracts.
Do not change runtime behavior.
Do not change Step 4 formulas.
Do not delete or modify legacy Step 3 behavior.
```

Recommended next behavior-change milestone after documentation:

```text
Generalize chromosome handling for 1-22, X, and Y.
```

## Final Principle

The pipeline contract is intentionally conservative.

The FAVOR-CLI refactor should modernize annotation and improve accessibility while preserving the current scientific prioritization and final GDS output contract.

Do not conflate annotation modernization with prioritization redesign.
