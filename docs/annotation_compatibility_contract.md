# IMPACT-SNV Annotation Compatibility Contract

## Purpose

This document defines the annotation compatibility contract required for IMPACT-SNV Step 4 prioritization.

It exists to support the staged refactor from the legacy FAVORannotator-based annotation step to a future FAVOR-CLI-backed annotation adapter.

The core purpose of this contract is to answer:

> What annotation fields must exist in an annotated GDS file before Step 4 can safely compute IMPACT prioritization scores and produce final `*_SNV_IMPACT.gds` outputs?

This contract is intentionally conservative. It describes the current minimum annotation interface expected by Step 4 and the compatibility requirements that any new FAVOR-CLI annotation path must satisfy.

This document should be read alongside:

- `.github/copilot-instructions.md`
- `docs/favorcli_refactor_plan.md`
- `docs/pipeline_contract.md`
- `docs/output_gds_contract.md`
- `docs/chromosome_handling_contract.md`

## Scope

This document covers:

- Required GDS annotation nodes consumed by Step 4.
- How each required annotation field is used in prioritization.
- Expected types and value conventions.
- Required fallbacks and failure behavior.
- Requirements for future FAVOR-CLI schema normalization.
- Requirements for GDS annotation injection.
- Validation requirements before Step 4 runs.

This document does **not** define:

- A new prioritization algorithm.
- New annotation-derived tiers.
- A replacement for `GeneList.txt`.
- A replacement for SeqArray GDS.
- A full FAVOR-CLI schema mapping.

The FAVOR-CLI output schema must be discovered empirically before final adapter mappings are implemented.

## Current Step 4 Annotation Dependency Summary

Step 4 currently expects an annotated SeqArray GDS file containing legacy FAVOR-style annotation fields under:

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

These fields form the initial **IMPACT annotation compatibility schema**.

Any FAVOR-CLI-backed annotation path must either:

1. Write these fields directly into the GDS under the same paths, or
2. Write a validated compatibility equivalent that Step 4 can consume without changing the scoring formulas.

The preferred initial strategy is to preserve the current paths so Step 4 can remain mostly unchanged.

## Required Annotation Namespace

Required compatibility namespace:

```text
annotation/info/FunctionalAnnotation/
```

The FAVOR-CLI adapter should inject normalized compatibility fields into this namespace so that current Step 4 code can continue to use `seqGetData()` against the expected paths.

Future implementations may additionally store richer FAVOR-CLI-native annotation fields elsewhere, but that should not replace the compatibility namespace until Step 4 is intentionally redesigned.

## Required Field Contract

## 1. `clnsig`

### Current GDS path

```text
annotation/info/FunctionalAnnotation/clnsig
```

### Purpose

`clnsig` stores ClinVar clinical significance labels.

It is used for:

- Tier 1 scoring.
- ClinVar pathogenic / likely pathogenic evidence.
- Final ClinVar boolean flag generation during post-processing.

### Tier dependency

Required for:

```text
Tier 1
ClinVar flag output
```

### Expected type

Character-like vector.

Expected length:

```text
length(clnsig) == number of variants in GDS
```

### Current expected value conventions

Step 4 currently depends on detecting pathogenic or likely pathogenic ClinVar values.

Relevant examples include:

```text
Pathogenic
Likely_pathogenic
Likely pathogenic
```

The legacy code has historically used exact or semi-exact token matching, so the adapter should normalize ClinVar labels conservatively.

### Required normalization behavior

The adapter should normalize ClinVar values so that Tier 1 behavior is preserved.

Recommended normalization rules:

```text
Pathogenic                → Pathogenic
pathogenic                → Pathogenic
Likely pathogenic         → Likely_pathogenic
Likely_pathogenic         → Likely_pathogenic
likely pathogenic         → Likely_pathogenic
likely_pathogenic         → Likely_pathogenic
```

For multi-valued ClinVar strings, the adapter should preserve enough delimiter information for Step 4 or post-processing to detect individual categories.

Supported delimiters should include, where practical:

```text
,
;
/
|
_
```

However, delimiter normalization must be tested against existing Step 4 and `seqarray_append.R` behavior.

### Missing value behavior

If ClinVar information is unavailable, the adapter should write an explicit fallback such as:

```text
not provided
```

or:

```text
NA
```

The selected fallback must be documented and consistently handled.

Preferred fallback:

```text
not provided
```

### Failure behavior

Missing `clnsig` should not necessarily block structural GDS validity, but it affects Tier 1 sensitivity and ClinVar flag output.

Recommended behavior:

- If the source ClinVar field is missing entirely, create `clnsig` with explicit fallback values.
- Record fallback usage in adapter output reports and GDS provenance.
- Do not silently omit the node.

### Scientific risk

Incorrect ClinVar normalization can change Tier 1 counts.

Any change to Tier 1 ClinVar parsing should be tested with known pathogenic, likely pathogenic, benign, uncertain significance, and missing examples.

## 2. `genecode_comprehensive_info`

### Current GDS path

```text
annotation/info/FunctionalAnnotation/genecode_comprehensive_info
```

### Purpose

`genecode_comprehensive_info` stores gene or transcript annotation information used to extract gene symbols.

Step 4 uses this field to match variants against `GeneList.txt`.

### Tier dependency

Required for all scoring tiers because all tiers depend on a non-zero gene-disease association score from `GeneList.txt`.

Required for:

```text
Tier 1
Tier 2
Tier 3
Tier 4
```

### Expected type

Character-like vector.

Expected length:

```text
length(genecode_comprehensive_info) == number of variants in GDS
```

### Expected value conventions

The current Step 4 logic expects gene symbols to be extractable from this field.

Examples may include:

```text
GJB2
OTOF
MYO6
NONE
NONE(dist=NONE)
GENE1,GENE2
GENE1(dist=123)
```

The current extraction logic is designed around legacy FAVORannotator-style gene annotation strings.

### Required normalization behavior

The FAVOR-CLI adapter should preserve gene-symbol information in a form compatible with current Step 4 extraction logic.

Preferred compatibility behavior:

- Preserve simple HGNC gene symbols where available.
- Preserve comma-separated multi-gene entries where applicable.
- Preserve or synthesize legacy-compatible `dist=` formatting only if required by current parsing behavior.
- Avoid inventing gene symbols.

If FAVOR-CLI provides a simpler gene symbol field, the adapter may map it into this field if the semantics are clear and documented.

### Missing value behavior

If gene annotation is unavailable, use:

```text
NONE
```

### Failure behavior

If gene information is absent for all variants, Step 4 may produce zero scores for all variants.

Recommended behavior:

- Write `NONE` fallback values if the source field is missing.
- Record fallback usage.
- Warn clearly that gene matching and scoring may be disabled or severely reduced.

### Scientific risk

Incorrect gene-symbol mapping can assign phenotype relevance to the wrong variant or fail to score a relevant variant.

This is one of the highest-risk fields in the adapter.

## 3. `genecode_comprehensive_exonic_category`

### Current GDS path

```text
annotation/info/FunctionalAnnotation/genecode_comprehensive_exonic_category
```

### Purpose

Stores GENCODE-based exonic consequence categories.

Used to identify variants eligible for Tier 2 and Tier 3 scoring.

### Tier dependency

Required for:

```text
Tier 2
Tier 3
```

### Expected type

Character-like vector.

Expected length:

```text
length(genecode_comprehensive_exonic_category) == number of variants in GDS
```

### Current expected value conventions

Tier 2 consequence labels include:

```text
frameshift insertion
frameshift deletion
stopgain
```

Tier 3 consequence labels include:

```text
nonsynonymous SNV
nonframeshift deletion
nonframeshift insertion
stoploss
```

### Required normalization behavior

If FAVOR-CLI uses Sequence Ontology-like labels, the adapter may normalize them to legacy labels only when the mapping is scientifically valid.

Possible mappings to validate:

```text
stop_gained          → stopgain
frameshift_variant   → frameshift insertion OR frameshift deletion, if indel direction is known
missense_variant     → nonsynonymous SNV
inframe_deletion     → nonframeshift deletion
inframe_insertion    → nonframeshift insertion
stop_lost            → stoploss
```

Ambiguous mappings must fail validation or be explicitly reported.

For example:

```text
frameshift_variant
```

is not always enough by itself to distinguish insertion from deletion unless REF/ALT allele length is used.

### Missing value behavior

If unavailable, use an explicit missing or neutral value such as:

```text
NA
```

or:

```text
NONE
```

Preferred fallback:

```text
NONE
```

### Failure behavior

If this field is missing, Step 4 Tier 2/3 scoring can still potentially use RefSeq or UCSC exonic categories if present.

Recommended behavior:

- Create the node with fallback values if the source field is unavailable.
- Record fallback usage.
- Do not silently copy another transcript system into this field unless explicitly documented.

### Scientific risk

Copying a single consequence source into GENCODE, RefSeq, and UCSC fields can overstate evidence agreement across transcript models.

If the adapter must do this temporarily, it must be documented as a compatibility fallback and recorded in provenance.

## 4. `refseq_exonic_category`

### Current GDS path

```text
annotation/info/FunctionalAnnotation/refseq_exonic_category
```

### Purpose

Stores RefSeq-based exonic consequence categories.

Used to identify variants eligible for Tier 2 and Tier 3 scoring.

### Tier dependency

Required for:

```text
Tier 2
Tier 3
```

### Expected type

Character-like vector.

Expected length:

```text
length(refseq_exonic_category) == number of variants in GDS
```

### Expected value conventions

Same Tier 2 and Tier 3 consequence labels as above:

Tier 2:

```text
frameshift insertion
frameshift deletion
stopgain
```

Tier 3:

```text
nonsynonymous SNV
nonframeshift deletion
nonframeshift insertion
stoploss
```

### Required normalization behavior

The adapter should map FAVOR-CLI RefSeq-specific consequence fields into this node only if a RefSeq-specific source field exists or the mapping is explicitly documented as a fallback.

### Missing value behavior

Preferred fallback:

```text
NONE
```

### Failure behavior

The node should still exist even if no RefSeq-specific annotation is available.

Recommended behavior:

- Fill with `NONE`.
- Record fallback usage.
- Do not silently omit the node.

### Scientific risk

If RefSeq consequences are unavailable and GENCODE consequences are copied into this field, the resulting Tier 2/3 classification may appear more strongly supported than it is.

Any such fallback must be clearly flagged.

## 5. `ucsc_exonic_category`

### Current GDS path

```text
annotation/info/FunctionalAnnotation/ucsc_exonic_category
```

### Purpose

Stores UCSC-based exonic consequence categories.

Used to identify variants eligible for Tier 2 and Tier 3 scoring.

### Tier dependency

Required for:

```text
Tier 2
Tier 3
```

### Expected type

Character-like vector.

Expected length:

```text
length(ucsc_exonic_category) == number of variants in GDS
```

### Expected value conventions

Same Tier 2 and Tier 3 consequence labels as above.

Tier 2:

```text
frameshift insertion
frameshift deletion
stopgain
```

Tier 3:

```text
nonsynonymous SNV
nonframeshift deletion
nonframeshift insertion
stoploss
```

### Required normalization behavior

The adapter should map FAVOR-CLI UCSC-specific consequence fields into this node only if a UCSC-specific source field exists or the mapping is explicitly documented as a fallback.

### Missing value behavior

Preferred fallback:

```text
NONE
```

### Failure behavior

The node should still exist even if no UCSC-specific annotation is available.

Recommended behavior:

- Fill with `NONE`.
- Record fallback usage.
- Do not silently omit the node.

### Scientific risk

Same as RefSeq: copying consequence values from another transcript model into this node can overstate evidence.

## 6. `apc_protein_function_v3`

### Current GDS path

```text
annotation/info/FunctionalAnnotation/apc_protein_function_v3
```

### Purpose

Stores APC protein function score used for Tier 4 prioritization.

### Tier dependency

Required for:

```text
Tier 4
```

### Expected type

Numeric vector.

Expected length:

```text
length(apc_protein_function_v3) == number of variants in GDS
```

### Expected value conventions

The current Step 4 logic treats missing values as zero and caps APC protein function values before normalization.

Current conceptual flow:

```text
NA → 0
cap values at 40
normalized_APC = APC / 40
```

### Required normalization behavior

The FAVOR-CLI adapter must verify that any source APC field is semantically equivalent to legacy `apc_protein_function_v3` before mapping it.

Potential source candidates may include names similar to:

```text
apc_protein_function_v3
aPC_protein_function_v3
apc_protein_function
```

However, field names alone are not sufficient. The scale and definition must be compatible.

### Missing value behavior

If APC protein function is unavailable, use:

```text
0
```

and record fallback usage.

### Failure behavior

If an APC-like source field exists but its scale or meaning is unclear, the adapter should fail fast or require explicit configuration rather than silently mapping it.

### Scientific risk

Incorrect APC score mapping can change Tier 4 prioritization.

This is a high-risk field and requires empirical validation against real FAVOR-CLI output and documentation.

## Gene-Disease Association Contract

Step 4 also requires a phenotype-specific gene-disease association file.

Current expected file:

```text
GeneList.txt
```

Required columns:

```text
symbol
globalScore
```

This contract must not change during the FAVOR-CLI adapter refactor.

The annotation adapter must ensure that gene symbols extracted from `genecode_comprehensive_info` can match the `symbol` values in `GeneList.txt`.

## Variant Identity and Row Alignment Contract

Annotation values must align exactly with the variant order in the GDS.

This is a critical scientific correctness requirement.

The adapter should maintain a variant identity table with at least:

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

Before writing annotations into GDS, validate:

```text
number of normalized annotation rows == number of GDS variants
all canonical keys are unique or duplicates are explicitly resolved
row order matches GDS variant order
chromosome naming is normalized
REF/ALT alleles match expected GDS alleles
```

If row alignment cannot be proven, annotation injection must fail.

## Chromosome Compatibility

The annotation compatibility schema must work for:

```text
1-22, X, Y
```

Chromosome aliases should be normalized explicitly:

```text
chr1 → 1
chrX → X
chrY → Y
```

The adapter must not silently convert:

```text
23 → X
24 → Y
```

unless that behavior is explicitly configured and recorded in provenance.

Mitochondrial support is not required for the initial FAVOR-CLI adapter.

## FAVOR-CLI Schema Discovery Requirement

The exact FAVOR-CLI output schema must not be guessed.

Before implementing final mappings, collect real FAVOR-CLI output on a controlled small variant set and document:

```text
FAVOR-CLI version
FAVOR database version
command used
input variant format
output format
output column names
sample output rows
candidate mappings to required IMPACT fields
unknown or unsafe mappings
```

Each required compatibility field should be classified as one of:

```text
direct match
rename only
transformation required
fallback required
unavailable
unsafe / unknown
```

Only fields classified as direct match, rename only, or validated transformation should be used without manual review.

## Adapter Mapping Configuration

The FAVOR-CLI adapter should use a versioned mapping configuration file rather than hard-coding all source field names directly into R or Bash scripts.

Suggested file:

```text
step3_favorcli_annotation/src/schema_mapping.yml
```

Suggested mapping structure:

```yaml
schema_version: impact-snv-annotation-schema-1.0.0

fields:
  clnsig:
    required: true
    type: character
    fallback: not provided
    source_candidates:
      - clnsig
      - clinvar_clnsig
      - clinical_significance

  genecode_comprehensive_info:
    required: true
    type: character
    fallback: NONE
    source_candidates:
      - genecode_comprehensive_info
      - gencode_gene_symbol
      - gene_symbol

  apc_protein_function_v3:
    required: true
    type: numeric
    fallback: 0
    source_candidates:
      - apc_protein_function_v3
      - aPC_protein_function_v3
```

The actual source candidates must be updated after real FAVOR-CLI schema discovery.

## GDS Injection Contract

Normalized adapter output should be injected into the GDS under:

```text
annotation/info/FunctionalAnnotation/
```

Required output nodes:

```text
annotation/info/FunctionalAnnotation/clnsig
annotation/info/FunctionalAnnotation/genecode_comprehensive_info
annotation/info/FunctionalAnnotation/genecode_comprehensive_exonic_category
annotation/info/FunctionalAnnotation/refseq_exonic_category
annotation/info/FunctionalAnnotation/ucsc_exonic_category
annotation/info/FunctionalAnnotation/apc_protein_function_v3
```

After injection, each node must be readable with `seqGetData()`.

The GDS injection process should also store provenance under:

```text
annotation/info/IMPACT_AnnotationProvenance/
annotation/info/IMPACT_AnnotationCompatibility/
```

## Fallback Reporting Contract

Fallbacks must never be silent.

For each field, the adapter should report:

```text
field name
source field used
fallback value
number of fallback values written
percentage of variants using fallback
reason for fallback
whether fallback affects scoring
```

Fallback reports may be written as:

```text
fallback_summary.tsv
schema_validation_report.txt
annotation_adapter_report.json
```

The final format can be chosen during implementation, but fallback reporting is required.

## Pre-Step-4 Schema Validation

Before Step 4 runs, validate required nodes.

For each required node, check:

```text
exists
readable with seqGetData()
expected type
length equals GDS variant count
missingness/fallback status is available
```

For variant alignment, check:

```text
variant count matches
variant keys align
no unresolved duplicate keys
chromosome naming is normalized
```

If validation fails, Step 4 should not run.

## Compatibility with Legacy FAVORannotator GDS Files

Step 4 must continue to support legacy FAVORannotator-produced GDS files during the transition period.

The FAVOR-CLI adapter should target compatibility with the legacy annotation namespace rather than forcing immediate changes to Step 4.

Legacy compatibility means:

- Existing `annotation/info/FunctionalAnnotation/*` nodes remain valid Step 4 inputs.
- Step 4 formulas remain unchanged.
- Final output nodes remain unchanged.
- Final output naming remains unchanged.

## Scientific Correctness Risks

The highest-risk compatibility areas are:

1. Variant row-order alignment between FAVOR-CLI output and GDS variants.
2. Gene-symbol mapping into `genecode_comprehensive_info`.
3. ClinVar significance normalization for Tier 1.
4. Consequence normalization for Tier 2 and Tier 3.
5. APC protein function score equivalence for Tier 4.
6. Copying one transcript annotation source into multiple transcript-specific legacy fields.
7. Silent fallbacks that allow the pipeline to run but reduce scientific sensitivity.

These risks should be addressed with explicit validation and golden-dataset testing.

## Testing Expectations

Adapter and schema tests should include:

- Missing ClinVar field.
- Missing gene annotation field.
- Missing GENCODE consequence field.
- Missing RefSeq consequence field.
- Missing UCSC consequence field.
- Missing APC score.
- Unknown ClinVar labels.
- Unknown consequence labels.
- X chromosome variant.
- Y chromosome variant.
- SNV.
- Insertion.
- Deletion.
- Duplicate variant key.
- Row-order mismatch.
- Incorrect REF/ALT mapping.

Golden validation should include variants that exercise:

```text
Tier 1
Tier 2
Tier 3
Tier 4
No-score controls
```

## Implementation Guidance for Coding Agents

When implementing code against this contract:

- Do not modify Step 4 formulas.
- Do not remove legacy FAVORannotator support.
- Do not assume FAVOR-CLI schema before schema discovery.
- Implement validation before relying on new annotation files.
- Prefer explicit mapping configuration over hard-coded assumptions.
- Fail fast on unsafe mappings.
- Report all fallbacks.
- Preserve final output compatibility.

## Future Revisions

This contract may be revised after real FAVOR-CLI schema discovery.

Any revision should clearly document:

- The FAVOR-CLI version tested.
- The FAVOR database version tested.
- The output schema observed.
- Which fields map directly.
- Which fields require transformation.
- Which legacy fields cannot be safely reproduced.
- Whether Step 4 must be minimally adapted.

## Final Principle

The FAVOR-CLI adapter is allowed to modernize annotation sourcing.

It is not allowed to silently change the scientific interpretation contract of IMPACT-SNV.

If a mapping is uncertain, report it and fail safely rather than producing misleading scores.
