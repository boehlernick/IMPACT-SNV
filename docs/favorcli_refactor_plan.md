# IMPACT-SNV FAVOR-CLI Refactor Plan

## Purpose

This document defines the intended staged refactor of **IMPACT-SNV** to replace the legacy FAVORannotator-based annotation step with a future **FAVOR-CLI-backed annotation adapter**, while preserving the existing IMPACT-SNV prioritization logic and final GDS output contract.

This is a planning and implementation guide for human developers and GitHub Copilot coding agents. It should be read alongside:

- `.github/copilot-instructions.md`
- `docs/pipeline_contract.md` *(to be created)*
- `docs/annotation_compatibility_contract.md` *(to be created)*
- `docs/output_gds_contract.md` *(to be created)*

## Executive Summary

The recommended architecture is:

```text
normalized VCF/GDS variant identity
→ FAVOR-CLI annotation
→ normalized IMPACT annotation compatibility schema
→ inject annotations into SeqArray GDS
→ existing IMPACT prioritization
→ final *_SNV_IMPACT.gds
```

The refactor should treat **FAVOR-CLI as a replaceable annotation backend**, not as a replacement for the full IMPACT-SNV workflow.

The immediate goal is **not** to redesign IMPACT-SNV's scoring method. The immediate goal is to modernize annotation and improve accessibility while preserving the current downstream contract.

## Current Pipeline

The current repository is organized as a four-step workflow:

```text
VCF files
→ step1_vcf_merge/
→ step2_vcf2gds/
→ Step3_favorannotator-rap/
→ step4_impact_prioritization/
→ final {sample_id}_SNV_IMPACT.gds
```

### Step 1: VCF Merge and Chromosome Split

Current folder:

```text
step1_vcf_merge/
```

Current responsibilities:

- Accept one or more VCF or VCF.GZ files.
- Download inputs in the DNAnexus/RAP environment.
- Ensure VCF files are bgzipped and indexed.
- Keep only `GT` FORMAT fields.
- Merge input VCFs.
- Normalize and split multiallelic variants using `bcftools norm`.
- Split merged VCF by chromosome.

Current implementation note:

- Chromosome splitting now supports `1-22`, `X`, and `Y`.
- Existing autosome behavior remains the compatibility baseline that later milestones must preserve.

### Step 2: VCF to GDS Conversion

Current folder:

```text
step2_vcf2gds/
```

Current responsibilities:

- Convert chromosome-specific VCF files to SeqArray GDS files.
- Use `SeqArray::seqVCF2GDS`.
- Add `annotation/info/QC_label = PASS` to all variants.

This step should remain mostly stable during the initial FAVOR-CLI refactor.

### Step 3: Legacy FAVORannotator Annotation

Current folder:

```text
Step3_favorannotator-rap/
```

Current responsibilities:

- Accept chromosome-level GDS files.
- Use a legacy FAVORannotator-style R workflow.
- Join GDS variant identities against FAVOR CSV database chunks.
- Write FAVOR-style annotations into the GDS under:

```text
annotation/info/FunctionalAnnotation/
```

Known current limitations:

- Tightly coupled to DNAnexus/RAP execution.
- Tightly coupled to legacy FAVORannotator behavior.
- Uses hard-coded/static FAVOR database assumptions.
- Autosome-oriented assumptions are present.
- Does not provide a clean abstraction for swapping annotation backends.

### Step 4: IMPACT Variant Prioritization

Current folder:

```text
step4_impact_prioritization/
```

Current responsibilities:

- Read annotated chromosome-level GDS files.
- Read `GeneList.txt` with phenotype-specific gene-disease association scores.
- Compute IMPACT pathogenicity scores and tier assignments.
- Export per-sample chromosome-level GDS files.
- Merge per-sample chromosome files into final:

```text
{sample_id}_SNV_IMPACT.gds
```

- Post-process final GDS files to add:

```text
annotation/info/impact_score
annotation/info/impact_score_calc
annotation/info/tier
annotation/info/clnsig_flags/*
```

## Refactor Goals

The refactor should accomplish the following:

1. Replace legacy FAVORannotator with a FAVOR-CLI-backed annotation path.
2. Allow users to leverage current FAVOR database releases and explicit FAVOR database configuration.
3. Preserve input VCF expectations where practical.
4. Preserve the final per-sample `*_SNV_IMPACT.gds` output contract.
5. Preserve current IMPACT prioritization formulas.
6. Add X and Y chromosome annotation support.
7. Reduce RAP/DNAnexus-specific coupling where practical.
8. Keep legacy FAVORannotator support available until the FAVOR-CLI path is validated.
9. Add validation and provenance so annotation differences can be explained and audited.

## Non-Goals for the Initial Refactor

Do not include the following in the initial implementation unless explicitly requested:

- Redesigning IMPACT prioritization tiers.
- Replacing GDS with Parquet, VCF, or another downstream exchange format.
- Introducing tissue-specific FAVOR annotations into scoring.
- Adding ACMG-like interpretation logic.
- Replacing IMPACT-VIS compatibility requirements.
- Introducing a required workflow engine such as Nextflow, Snakemake, or WDL.
- Removing `Step3_favorannotator-rap/`.
- Changing `GeneList.txt` format.
- Changing final output naming.

## Recommended Architecture

Use an adapter-based architecture:

```text
GDS or normalized VCF variant identity
→ FAVOR-CLI raw annotation output
→ IMPACT annotation compatibility adapter
→ normalized annotation table
→ GDS injection
→ Step 4 prioritization
```

This architecture isolates annotation-tool changes from Step 4's scientific prioritization logic.

### Why an Adapter Is Preferred

The current Step 4 is tightly coupled to a small set of legacy annotation fields under:

```text
annotation/info/FunctionalAnnotation/
```

Rather than rewriting Step 4 immediately, the safer path is to make FAVOR-CLI output look like the legacy internal annotation contract.

This allows:

- Existing prioritization formulas to remain stable.
- Existing final GDS output structure to remain stable.
- Legacy and new annotation paths to coexist during validation.
- Annotation schema changes to be isolated in one adapter layer.

## Annotation Compatibility Contract

The adapter must provide, or explicitly validate the absence of, the fields Step 4 currently consumes.

Required compatibility fields:

```text
annotation/info/FunctionalAnnotation/clnsig
annotation/info/FunctionalAnnotation/genecode_comprehensive_info
annotation/info/FunctionalAnnotation/genecode_comprehensive_exonic_category
annotation/info/FunctionalAnnotation/refseq_exonic_category
annotation/info/FunctionalAnnotation/ucsc_exonic_category
annotation/info/FunctionalAnnotation/apc_protein_function_v3
```

### Field Responsibilities

#### `clnsig`

Used for:

- Tier 1 ClinVar pathogenic / likely pathogenic evidence.
- Final ClinVar flag post-processing.

Compatibility requirements:

- Must be character-like.
- Must preserve or normalize clinical significance labels into a form Step 4 can parse.
- Missing values should be explicit, for example `not provided`, rather than silently absent.

#### `genecode_comprehensive_info`

Used for:

- Extracting gene symbols.
- Matching variants to phenotype-associated genes in `GeneList.txt`.

Compatibility requirements:

- Must be character-like.
- Must preserve gene-symbol information in a format compatible with current extraction logic, or Step 4 must be minimally adapted with equivalent parsing.
- Missing values should use an explicit fallback such as `NONE`.

#### Exonic category fields

Required fields:

```text
annotation/info/FunctionalAnnotation/genecode_comprehensive_exonic_category
annotation/info/FunctionalAnnotation/refseq_exonic_category
annotation/info/FunctionalAnnotation/ucsc_exonic_category
```

Used for:

- Tier 2 consequence detection.
- Tier 3 consequence detection.

Compatibility requirements:

- Must be character-like.
- Must normalize consequence labels to current expected values where appropriate.
- Must not overstate evidence by pretending that one annotation source is equivalent to all three transcript systems unless explicitly documented.

Current expected Tier 2 labels include:

```text
frameshift insertion
frameshift deletion
stopgain
```

Current expected Tier 3 labels include:

```text
nonsynonymous SNV
nonframeshift deletion
nonframeshift insertion
stoploss
```

Possible normalization examples:

```text
stop_gained          → stopgain
missense_variant     → nonsynonymous SNV
inframe_deletion     → nonframeshift deletion
inframe_insertion    → nonframeshift insertion
stop_lost            → stoploss
```

These mappings must be validated before being treated as scientifically equivalent.

#### `apc_protein_function_v3`

Used for:

- Tier 4 APC protein-function evidence.

Compatibility requirements:

- Must be numeric.
- Missing values should fall back to `0` only with explicit reporting.
- The adapter must validate that the source score is semantically and numerically compatible with the legacy APC protein function score expected by Step 4.

## Variant Identity Contract

The adapter must preserve a stable join between GDS variants and FAVOR-CLI annotations.

Required variant identity fields should include:

```text
variant.id
chromosome
position
ref
alt
canonical_variant_key
```

The canonical variant key should be documented and consistently generated.

Recommended key format:

```text
chromosome-position-ref-alt
```

Examples:

```text
1-12345-A-G
X-54321-C-T
Y-22222-G-A
```

Chromosome normalization must be explicit.

Recommended internal chromosome representation:

```text
1, 2, ..., 22, X, Y
```

Accepted external aliases may include:

```text
chr1 → 1
chrX → X
chrY → Y
```

Do not silently convert:

```text
23 → X
24 → Y
```

unless explicitly configured.

## Provenance Requirements

The FAVOR-CLI adapter and GDS injection layer should store provenance under stable namespaces such as:

```text
annotation/info/IMPACT_AnnotationProvenance/
annotation/info/IMPACT_AnnotationCompatibility/
```

Minimum recommended provenance fields:

```text
favor_cli_version
favor_database_version
favor_database_release_date
reference_genome_build
reference_genome_fasta_name
reference_genome_fasta_checksum
annotation_command
annotation_date_utc
adapter_version
adapter_schema_version
input_vcf_name
input_vcf_checksum
input_gds_name
input_gds_checksum
annotation_output_format
annotation_output_checksum
chromosome_normalization_policy
variant_key_policy
field_mapping_policy
fallback_policy
```

The purpose of provenance is to distinguish:

- Expected annotation changes due to updated FAVOR data.
- Adapter or schema mapping errors.
- Pipeline regressions.

## Pre-Step-4 Validation Requirements

Before Step 4 runs on adapter-injected GDS files, validate that the required annotation contract is satisfied.

For each required node, validate:

```text
exists
readable by seqGetData()
expected type
length == GDS variant count
missingness/fallback status recorded
```

Also validate:

```text
variant row order is aligned
no duplicate unresolved variant keys
no dropped variants unless explicitly reported
chromosome labels are normalized
X/Y variants are retained when present
```

If validation fails, the pipeline should fail fast with a clear error message.

## Milestone Roadmap

The refactor should be implemented as small, reviewable PRs or milestones.

## Milestone 1: Documentation and Contracts

### Objective

Document the current pipeline contract, annotation compatibility contract, and final output contract before runtime behavior changes.

### Files to add or update

```text
.github/copilot-instructions.md
docs/favorcli_refactor_plan.md
docs/pipeline_contract.md
docs/annotation_compatibility_contract.md
docs/output_gds_contract.md
README.md
```

### Required content

Document:

- Current four-step workflow.
- Current Step 4 required annotation nodes.
- Current prioritization formulas.
- Final GDS output contract.
- Current limitations.
- Intended FAVOR-CLI adapter architecture.

### Acceptance criteria

- No runtime code changes.
- No scoring formula changes.
- No file deletion.
- Documentation clearly distinguishes current behavior from planned behavior.
- A reviewer can understand what must remain stable during the refactor.

## Milestone 2: Generalize Chromosome Handling

### Objective

Support chromosomes:

```text
1-22, X, Y
```

while preserving existing autosome behavior.

### Likely files

```text
step1_vcf_merge/src/IMPACT_SNV_indel.py
step1_vcf_merge/Readme.md
step4_impact_prioritization/resources/home/dnanexus/IMPACT-prioritization.r
step4_impact_prioritization/README.md
```

### Implementation requirements

- Replace autosome-only splitting with configurable chromosome handling.
- Preserve output naming conventions:

```text
merged_chr1.vcf.gz
merged_chrX.vcf.gz
merged_chrY.vcf.gz
merged_chr1.gds
merged_chrX.gds
merged_chrY.gds
```

- Ensure Step 4 can discover and merge per-sample X/Y files:

```text
{sample_id}_chrX.gds
{sample_id}_chrY.gds
```

- Do not introduce mitochondrial chromosome support as required behavior yet.

### Acceptance criteria

- Existing autosome-only runs still work.
- X/Y VCF split outputs are generated when input variants exist.
- Step 4 can process X/Y chromosome GDS files.
- Step 4 formulas are unchanged.

## Milestone 3: FAVOR-CLI Step 3 Skeleton

### Objective

Add a new experimental Step 3 module for FAVOR-CLI annotation without replacing legacy FAVORannotator yet.

### New folder

```text
step3_favorcli_annotation/
```

Suggested structure:

```text
step3_favorcli_annotation/
  README.md
  dxapp.json
  src/code.sh
  src/run_favorcli.sh
  src/extract_variant_identity.R
  src/validate_favorcli_config.sh
  src/validate_favorcli_output.R
  tests/
```

### Functional requirements

- Accept a GDS file as input.
- Extract canonical variant identity from GDS.
- Produce a FAVOR-CLI-ready input file or canonical intermediate TSV.
- Support dry-run mode.
- Make configuration explicit:

```text
FAVOR-CLI executable path
FAVOR database path
FAVOR database version
reference genome build
threads
memory budget
output directory
dry-run mode
```

### Acceptance criteria

- Legacy Step 3 remains untouched.
- Dry-run mode works without requiring FAVOR-CLI execution.
- Variant identity table is produced.
- No GDS mutation occurs in this milestone.

## Milestone 4: Annotation Schema Adapter

### Objective

Convert FAVOR-CLI raw output into the normalized IMPACT annotation compatibility schema required by Step 4.

### Likely files

```text
step3_favorcli_annotation/src/schema_mapping.yml
step3_favorcli_annotation/src/normalize_favorcli_annotations.R
step3_favorcli_annotation/src/validate_impact_annotation_schema.R
step3_favorcli_annotation/tests/
docs/annotation_compatibility_contract.md
```

### Implementation requirements

- Use a versioned mapping file.
- Normalize into IMPACT-owned compatibility fields.
- Validate all required fields.
- Validate row counts against the variant identity table.
- Report fallback usage explicitly.
- Do not assume real FAVOR-CLI schema unless real output examples are available.

Required normalized fields:

```text
clnsig
genecode_comprehensive_info
genecode_comprehensive_exonic_category
refseq_exonic_category
ucsc_exonic_category
apc_protein_function_v3
```

### Acceptance criteria

- Adapter can process a fixture into a normalized annotation table.
- Schema validator fails clearly if required fields are missing.
- Fallbacks are reported.
- No Step 4 formulas are changed.

## Milestone 5: GDS Annotation Injection

### Objective

Write normalized adapter outputs into the GDS namespace consumed by Step 4:

```text
annotation/info/FunctionalAnnotation/
```

### Likely files

```text
step3_favorcli_annotation/src/inject_annotations_gds.R
step3_favorcli_annotation/src/validate_gds_annotations.R
step3_favorcli_annotation/src/code.sh
step3_favorcli_annotation/README.md
```

### Implementation requirements

- Open GDS in write mode using SeqArray/gdsfmt.
- Validate row count and variant key alignment before writing.
- Write required compatibility nodes under `annotation/info/FunctionalAnnotation/`.
- Validate nodes after writing.
- Store provenance metadata.

### Acceptance criteria

- Injected GDS passes pre-Step-4 schema validation.
- Existing Step 4 can read required annotation nodes.
- Provenance metadata is present.
- Row-order mismatches fail fast.

## Milestone 6: Minimal Step 4 Robustness Updates

### Objective

Make Step 4 robust enough to consume both legacy FAVORannotator GDS files and adapter-injected GDS files.

### Likely files

```text
step4_impact_prioritization/resources/home/dnanexus/IMPACT-prioritization.r
step4_impact_prioritization/resources/home/dnanexus/seqarray_append.R
step4_impact_prioritization/README.md
```

### Allowed changes

- Add preflight schema validation.
- Improve missing-node error messages.
- Improve ClinVar parsing while preserving existing behavior.
- Ensure X/Y file discovery and merging works.

### Forbidden changes

- Do not change Tier 1-4 formulas.
- Do not change final output names.
- Do not rename final output nodes.
- Do not redesign prioritization logic.

### Acceptance criteria

- Legacy annotated GDS still works.
- Adapter-injected GDS works.
- X/Y chromosome files can be processed.
- Final `{sample_id}_SNV_IMPACT.gds` files are produced.

## Milestone 7: Golden-Dataset Validation

### Objective

Create a small validation suite proving the refactor preserves structure and prioritization behavior.

### Suggested folder

```text
tests/golden/
```

Suggested contents:

```text
tests/golden/README.md
tests/golden/run_golden_validation.sh
tests/golden/fixtures/
tests/golden/expected/
tests/golden/reports/
```

### Required scenarios

Include variants that exercise:

```text
Tier 1
Tier 2
Tier 3
Tier 4
No-score controls
Autosomal variants
X chromosome variants
Y chromosome variants
SNVs
indels
missing ClinVar
missing APC
missing gene match
```

### Validation checks

Golden validation should check:

- Required GDS nodes exist.
- Node lengths match variant count.
- Variant identity is preserved.
- Genotypes are preserved for retained variants.
- X/Y variants are not dropped.
- Scores match expected formulas.
- Tiers are correct.
- Final `{sample_id}_SNV_IMPACT.gds` files are produced.

### Acceptance criteria

- A documented command runs the golden validation.
- Reports distinguish structural compatibility, scoring equivalence, and annotation-content drift.
- X/Y variants appear in final outputs when expected.

## Milestone 8: Portability Improvements

### Objective

Reduce RAP/DNAnexus-specific coupling while preserving RAP support.

### Possible additions

```text
scripts/run_step1_local.sh
scripts/run_step2_local.sh
scripts/run_step3_favorcli_local.sh
scripts/run_step4_local.sh
containers/Dockerfile.impact-snv
containers/Dockerfile.favorcli
workflow/local/
workflow/rap/
workflow/container/
```

### Implementation principles

- Separate platform wrappers from core logic.
- Keep DNAnexus applets as supported wrappers.
- Make local execution possible without `dx`.
- Do not bundle large FAVOR databases into containers by default.
- Treat FAVOR database path as an external mount/configuration.

### Acceptance criteria

- Core FAVOR-CLI adapter can run outside DNAnexus.
- RAP applet wrappers delegate to shared core scripts where practical.
- Documentation includes at least one local or containerized example.

## Milestone 9: Documentation and Migration Guide

### Objective

Update user-facing documentation for the refactored pipeline.

### Required documentation

- Supported input VCF assumptions.
- Reference genome requirements.
- Chromosome handling.
- FAVOR-CLI setup.
- FAVOR database version pinning.
- Adapter schema.
- Provenance metadata.
- Local/container/RAP execution examples.
- Migration from legacy FAVORannotator.
- Expected differences due to updated FAVOR annotations.

### Acceptance criteria

- README reflects the new architecture.
- Legacy FAVORannotator instructions remain available but clearly marked legacy.
- Users can understand how to run the new FAVOR-CLI path.
- Users can understand why updated annotations may change results.

## Milestone 10: Legacy FAVORannotator Deprecation Gate

### Objective

Only after validation, mark legacy FAVORannotator as deprecated.

### Requirements

- Do not delete legacy Step 3 immediately.
- Add a deprecation notice.
- Keep a version tag or branch where legacy behavior is preserved.
- Confirm Step 4 still accepts legacy annotated GDS files.
- Provide migration instructions.

### Acceptance criteria

- FAVOR-CLI adapter path is the default documented workflow.
- Legacy path remains available for users who need it.
- Golden validation supports the deprecation decision.

## First Coding Agent Task

The first actual coding-agent task should be limited to documentation and contracts.

Recommended first prompt:

```text
You are working in the IMPACT-SNV repository.

I am beginning a staged refactor to replace the legacy FAVORannotator-based annotation step with a future FAVOR-CLI-backed annotation adapter. This first task is documentation and repo setup only.

Do not change runtime behavior in this task.

Tasks:
1. Inspect the repository structure and confirm the current four-step workflow.
2. Add or update:
   - .github/copilot-instructions.md
   - docs/favorcli_refactor_plan.md
   - docs/pipeline_contract.md
   - docs/annotation_compatibility_contract.md
   - docs/output_gds_contract.md
3. Update the root README.md with links to these docs.
4. Document the current Step 4 required annotation nodes.
5. Document the final output contract.
6. Document the current prioritization formulas and state they must not change during the FAVOR-CLI adapter refactor.
7. Document known limitations and the intended adapter-based refactor direction.

Acceptance criteria:
- No runtime code changes.
- No scoring formula changes.
- No file deletion.
- Documentation clearly distinguishes current behavior from planned FAVOR-CLI adapter behavior.
- The PR is small and reviewable.
```

## Review Checklist for Each PR

Before merging each PR, confirm:

- Does this PR change runtime behavior?
- If yes, is that behavior change expected and documented?
- Does this PR alter Step 4 formulas?
- Does this PR alter final output naming?
- Does this PR alter final output node locations?
- Does this PR preserve legacy FAVORannotator compatibility?
- Does this PR add validation before relying on new assumptions?
- Are fallbacks explicit and reported?
- Are chromosome naming assumptions documented?
- Are X/Y variants preserved where expected?
- Are provenance fields added when annotations are generated or injected?

## Open Unknowns

The following must be resolved empirically before full implementation:

1. Exact FAVOR-CLI input format to use for this workflow.
2. Exact FAVOR-CLI output schema for required fields.
3. Whether FAVOR-CLI provides direct equivalents for:
   - `clnsig`
   - `genecode_comprehensive_info`
   - GENCODE exonic category
   - RefSeq exonic category
   - UCSC exonic category
   - `apc_protein_function_v3`
4. Whether newer FAVOR annotations preserve legacy value conventions.
5. Whether APC protein function score scale/version is equivalent.
6. Whether X/Y annotation behavior differs from autosomes.
7. Whether any RAP-specific constraints complicate FAVOR database mounting or caching.

## Recommended Immediate Developer Workflow

1. Commit `.github/copilot-instructions.md`.
2. Add this `docs/favorcli_refactor_plan.md` file.
3. Add the three contract documents.
4. Run the coding agent only on Milestone 1.
5. Review and merge documentation changes.
6. Start a separate branch for chromosome handling.
7. Generate a tiny test VCF and run FAVOR-CLI manually to capture real output.
8. Use that real FAVOR-CLI output to design the adapter mapping.

## Final Principle

This refactor is about replacing a brittle legacy annotation backend while preserving the scientific prioritization and final GDS output contract.

Do not conflate annotation modernization with prioritization redesign.
