# IMPACT-SNV FAVOR-CLI Schema Discovery Protocol

## Purpose

This document defines the protocol for discovering and documenting the actual FAVOR-CLI output schema before implementing the IMPACT-SNV FAVOR-CLI annotation adapter.

The purpose of schema discovery is to prevent developers or coding agents from guessing how FAVOR-CLI output maps to the legacy IMPACT-SNV annotation fields currently consumed by Step 4.

The central question is:

> Which FAVOR-CLI output fields can be safely mapped, transformed, or adapted into the legacy `annotation/info/FunctionalAnnotation/*` fields required by IMPACT-SNV Step 4?

This protocol should be completed before implementing production adapter mappings.

This document should be read alongside:

- `.github/copilot-instructions.md`
- `docs/favorcli_refactor_plan.md`
- `docs/pipeline_contract.md`
- `docs/annotation_compatibility_contract.md`
- `docs/output_gds_contract.md`
- `docs/chromosome_handling_contract.md`

## Scope

This document covers:

- How to run FAVOR-CLI on a tiny controlled variant set.
- What metadata to capture.
- What output files to preserve.
- How to inspect FAVOR-CLI output columns.
- How to compare FAVOR-CLI output against IMPACT-SNV's required annotation compatibility schema.
- How to classify candidate field mappings.
- What must be resolved before implementing the annotation adapter.

This document does **not** define:

- Final FAVOR-CLI adapter implementation.
- Final schema mapping values.
- GDS annotation injection logic.
- Step 4 scoring changes.
- A replacement for the annotation compatibility contract.

## Why Schema Discovery Is Required

The existing IMPACT-SNV Step 4 prioritization code expects legacy FAVORannotator-style annotations under:

```text
annotation/info/FunctionalAnnotation/
```

The minimum required fields are:

```text
annotation/info/FunctionalAnnotation/clnsig
annotation/info/FunctionalAnnotation/genecode_comprehensive_info
annotation/info/FunctionalAnnotation/genecode_comprehensive_exonic_category
annotation/info/FunctionalAnnotation/refseq_exonic_category
annotation/info/FunctionalAnnotation/ucsc_exonic_category
annotation/info/FunctionalAnnotation/apc_protein_function_v3
```

FAVOR-CLI may expose equivalent, renamed, transformed, richer, or incompatible fields. The adapter must be based on observed output and documented mappings rather than assumptions.

Schema discovery is especially important for:

- ClinVar significance labels.
- Gene symbol and transcript annotation fields.
- GENCODE, RefSeq, and UCSC exonic consequence categories.
- APC protein function scores.
- Chromosome naming.
- Variant identity keys.
- X and Y chromosome support.

## Required Discovery Output Directory

Schema discovery outputs should be stored under:

```text
tests/favorcli_schema_discovery/
```

Recommended structure:

```text
tests/favorcli_schema_discovery/
  README.md
  input/
    input_variants.tsv
    input_variants.vcf
    input_variants.vcf.gz
  commands/
    favorcli_version.txt
    favorcli_manifest.txt
    favorcli_schema_command.txt
    favorcli_ingest_command.txt
    favorcli_annotate_command.txt
  output/
    favorcli_raw_output.*
    favorcli_columns.txt
    favorcli_schema_output.*
    favorcli_manifest_output.*
  notes/
    schema_mapping_notes.md
    required_field_mapping_matrix.tsv
    unresolved_fields.md
    unsafe_mappings.md
```

The exact output file extensions may vary depending on FAVOR-CLI output format, for example:

```text
.parquet
.tsv
.csv
.json
```

## Required Metadata to Capture

Each schema discovery run must capture:

```text
FAVOR-CLI version
FAVOR database version
FAVOR database path or configuration
reference genome build
input file name
input file checksum
command used for ingest, if applicable
command used for annotation
output format
output file name
output file checksum
run date/time
local environment or container information
```

Recommended metadata file:

```text
tests/favorcli_schema_discovery/README.md
```

or machine-readable sidecar:

```text
tests/favorcli_schema_discovery/output/discovery_metadata.json
```

## Recommended Tiny Input Variant Set

The schema discovery input should be small but biologically representative.

Include variants that can potentially exercise the fields needed by IMPACT-SNV:

```text
Autosomal SNV
Autosomal insertion
Autosomal deletion
X chromosome SNV
Y chromosome SNV
Known or likely ClinVar pathogenic example if available
Variant with expected gene annotation
Variant expected to have no ClinVar annotation
Variant expected to have protein consequence annotation
```

The input should include at least:

```text
chr1 or 1
chrX or X
chrY or Y
```

If real variants with known expected annotations are unavailable, use the smallest controlled set possible and clearly mark expected annotations as unknown.

## Preferred Input Formats to Test

If feasible, test FAVOR-CLI with the same general input representation that IMPACT-SNV can generate.

Candidate inputs:

```text
VCF or VCF.GZ
canonical variant TSV
FAVOR-CLI ingested intermediate format
```

The current IMPACT-SNV pipeline naturally produces:

```text
VCF.GZ files
SeqArray GDS files
GDS-derived variant identity tables
```

For the adapter path, the most important intermediate is a GDS-derived variant identity table containing:

```text
variant.id
chromosome
position
ref
alt
canonical_variant_key
```

Recommended canonical key:

```text
chromosome-position-ref-alt
```

Examples:

```text
1-12345-A-G
X-54321-C-T
Y-22222-G-A
```

## Commands to Capture

Capture the exact command used to inspect FAVOR-CLI version:

```bash
favor --version
```

Capture any available schema or manifest command output if supported by the installed FAVOR-CLI version:

```bash
favor schema
favor manifest
```

Capture the exact ingest command, if using an ingest step:

```bash
favor ingest input_variants.vcf.gz
```

Capture the exact annotation command:

```bash
favor annotate input.ingested.parquet
```

These command examples are conceptual. The exact commands should reflect the installed FAVOR-CLI version and the selected workflow.

Store commands in:

```text
tests/favorcli_schema_discovery/commands/
```

## Output Column Inspection

After running FAVOR-CLI, create a file listing output columns:

```text
tests/favorcli_schema_discovery/output/favorcli_columns.txt
```

For TSV/CSV output, this can be produced by reading the header.

For Parquet output, use an appropriate local tool or script to list schema fields.

The column list should preserve:

```text
column name
column type if available
nullable/missingness if available
example value if safe and useful
```

## Required Field Mapping Matrix

Create a mapping matrix at:

```text
tests/favorcli_schema_discovery/notes/required_field_mapping_matrix.tsv
```

The matrix should include one row for each required IMPACT-SNV compatibility field:

```text
clnsig
genecode_comprehensive_info
genecode_comprehensive_exonic_category
refseq_exonic_category
ucsc_exonic_category
apc_protein_function_v3
```

Recommended columns:

```text
impact_field
required_gds_path
required_for
expected_type
candidate_favorcli_field
candidate_favorcli_type
mapping_classification
transformation_required
fallback_if_missing
scientific_risk
status
notes
```

## Mapping Classification Categories

Each required field must be classified as one of:

```text
direct_match
rename_only
transformation_required
fallback_required
unavailable
unsafe_unknown
```

### `direct_match`

Use when FAVOR-CLI provides a field with the same name and compatible semantics.

Example:

```text
FAVOR-CLI field: clnsig
IMPACT field: clnsig
```

This still requires type and value validation.

### `rename_only`

Use when FAVOR-CLI provides a semantically equivalent field with a different name.

Example:

```text
FAVOR-CLI field: clinical_significance
IMPACT field: clnsig
```

This requires documentation but may not require value transformation.

### `transformation_required`

Use when FAVOR-CLI provides relevant information, but values must be transformed.

Example:

```text
stop_gained → stopgain
missense_variant → nonsynonymous SNV
Likely pathogenic → Likely_pathogenic
```

Transformations require tests.

### `fallback_required`

Use when FAVOR-CLI does not provide the field or no safe mapping exists, but IMPACT-SNV can safely write a neutral fallback value.

Examples:

```text
missing clnsig → not provided
missing gene info → NONE
missing APC score → 0
```

Fallbacks must be reported.

### `unavailable`

Use when FAVOR-CLI does not provide a required concept and no safe fallback should be assumed.

This may block adapter implementation or require Step 4 changes.

### `unsafe_unknown`

Use when a field appears related but semantic equivalence is unclear.

Examples:

```text
An APC-like score with unknown scale
A generic consequence field mapped into all transcript-specific legacy fields
A gene symbol field that may represent nearest gene rather than affected gene
```

Unsafe mappings must not be used silently.

## Required Field-Specific Questions

## ClinVar Significance

Questions to answer:

```text
Does FAVOR-CLI provide ClinVar significance?
What is the field name?
What delimiter is used for multiple labels?
Are labels space-separated, underscore-separated, lowercase, or normalized?
Are missing values explicit?
Can Pathogenic and Likely pathogenic be detected reliably?
```

Required IMPACT field:

```text
annotation/info/FunctionalAnnotation/clnsig
```

## Gene / Transcript Annotation

Questions to answer:

```text
Does FAVOR-CLI provide gene symbols?
Does it provide GENCODE comprehensive info or an equivalent?
Does the field identify affected genes, nearest genes, transcript genes, or regulatory target genes?
Can values be converted into the legacy gene-string format expected by Step 4?
Can multiple genes be represented?
```

Required IMPACT field:

```text
annotation/info/FunctionalAnnotation/genecode_comprehensive_info
```

## Exonic Consequence Categories

Questions to answer:

```text
Does FAVOR-CLI provide exonic consequence categories?
Are categories transcript-source-specific?
Are GENCODE, RefSeq, and UCSC annotations separate?
Are Sequence Ontology terms used?
Can frameshift insertion and frameshift deletion be distinguished?
Can inframe insertion and inframe deletion be distinguished?
```

Required IMPACT fields:

```text
annotation/info/FunctionalAnnotation/genecode_comprehensive_exonic_category
annotation/info/FunctionalAnnotation/refseq_exonic_category
annotation/info/FunctionalAnnotation/ucsc_exonic_category
```

## APC Protein Function Score

Questions to answer:

```text
Does FAVOR-CLI provide apc_protein_function_v3?
Does it provide another APC or aPC protein-function field?
What is the scale?
Is the score capped or uncapped?
Is it comparable to the legacy score expected by Step 4?
Are missing values represented as NA, null, blank, or zero?
```

Required IMPACT field:

```text
annotation/info/FunctionalAnnotation/apc_protein_function_v3
```

## Variant Identity and Join Keys

Questions to answer:

```text
Does FAVOR-CLI preserve input row order?
Does FAVOR-CLI emit a variant key?
Does FAVOR-CLI emit chromosome, position, ref, and alt separately?
Does chromosome naming include chr prefixes?
How are multiallelic variants represented?
Are indels normalized consistently with IMPACT-SNV Step 1?
Are duplicate variant keys possible?
```

Required adapter behavior:

```text
Annotations must align exactly to GDS variant order before injection.
```

If row order is not preserved, the adapter must perform an explicit join using validated variant keys.

## X/Y Chromosome Schema Discovery

Schema discovery must include X and Y variants where possible.

Questions to answer:

```text
Does FAVOR-CLI annotate chrX / X variants?
Does FAVOR-CLI annotate chrY / Y variants?
Does FAVOR-CLI require chr prefixes?
Does output use chrX or X?
Are X/Y annotations returned in the same schema as autosomes?
Are any fields missing for X/Y compared with autosomes?
```

If X/Y behavior differs from autosomes, document it in:

```text
tests/favorcli_schema_discovery/notes/schema_mapping_notes.md
```

## Required Notes Files

## `schema_mapping_notes.md`

This file should summarize:

- FAVOR-CLI version tested.
- FAVOR database version tested.
- Input variant set used.
- Output format produced.
- High-level schema observations.
- Candidate field mappings.
- Known uncertainties.
- Recommended next implementation step.

## `unresolved_fields.md`

This file should list required IMPACT fields that could not yet be mapped.

For each unresolved field, include:

```text
required IMPACT field
why unresolved
candidate FAVOR-CLI fields considered
what evidence is missing
whether fallback is acceptable
recommended next action
```

## `unsafe_mappings.md`

This file should list fields that appear related but should not be used without further validation.

Examples:

```text
APC-like score with unknown scale
Generic gene field with unclear semantics
Generic consequence field copied into transcript-specific nodes
ClinVar field with incompatible labels
```

## Adapter Implementation Gate

Do not implement production adapter mappings until schema discovery has answered the following:

1. Which FAVOR-CLI fields map to the six required IMPACT compatibility fields?
2. Which fields require transformation?
3. Which fields require fallbacks?
4. Which mappings are unsafe or unresolved?
5. How variant row alignment will be guaranteed?
6. Whether X/Y variants are annotated with the same schema?
7. Whether APC score semantics are compatible?

If these questions cannot be answered, implement only mock/fixture-based adapter scaffolding and mark production mapping as blocked.

## Minimum Acceptance Criteria for Schema Discovery

Schema discovery is complete enough to begin adapter implementation when:

- FAVOR-CLI version is recorded.
- FAVOR database version or configuration is recorded.
- Exact commands are recorded.
- Raw FAVOR-CLI output is preserved.
- Output columns are listed.
- Required field mapping matrix is created.
- Each required IMPACT field is classified.
- Unsafe mappings are documented.
- X/Y behavior is tested or explicitly marked untested.
- Variant identity and row-order behavior are understood.

## Recommended Coding Agent Prompt for Schema Discovery

Use this prompt only after the initial documentation contracts are committed:

```text
Using docs/favorcli_schema_discovery.md, create the schema discovery scaffold only.

Do not implement production FAVOR-CLI adapter mappings yet.

Tasks:
1. Create tests/favorcli_schema_discovery/ with the documented folder structure.
2. Add README.md explaining how to run schema discovery manually.
3. Add placeholder files or templates for:
   - commands/favorcli_version.txt
   - commands/favorcli_ingest_command.txt
   - commands/favorcli_annotate_command.txt
   - output/favorcli_columns.txt
   - notes/schema_mapping_notes.md
   - notes/required_field_mapping_matrix.tsv
   - notes/unresolved_fields.md
   - notes/unsafe_mappings.md
4. Do not assume FAVOR-CLI output schema.
5. Do not modify runtime pipeline behavior.
```

## Final Principle

Do not guess the FAVOR-CLI schema.

Capture real output, classify mappings, document uncertainty, and only then implement the adapter.
