# IMPACT-SNV Output GDS Contract

## Purpose

This document defines the expected final output contract for IMPACT-SNV per-sample GDS files.

It exists to support the staged refactor from the legacy FAVORannotator-based annotation step to a future FAVOR-CLI-backed annotation adapter while preserving downstream compatibility with IMPACT-SNV and IMPACT-VIS expectations.

The core purpose of this contract is to answer:

> What must the final `{sample_id}_SNV_IMPACT.gds` files contain, and what must not change during the FAVOR-CLI refactor?

This contract is intentionally conservative. It describes the current final output expectations and the minimum structural compatibility requirements that future refactor work must preserve.

This document should be read alongside:

- `.github/copilot-instructions.md`
- `docs/favorcli_refactor_plan.md`
- `docs/pipeline_contract.md`
- `docs/annotation_compatibility_contract.md`
- `docs/chromosome_handling_contract.md`

## Scope

This document covers:

- Final per-sample GDS file naming.
- Required final GDS nodes.
- Required IMPACT score annotations.
- Required tier annotations.
- Required ClinVar flag annotations.
- Core SeqArray node expectations.
- Sample and genotype preservation requirements.
- Compatibility expectations for IMPACT-VIS.
- Validation expectations for final GDS outputs.

This document does **not** define:

- The FAVOR-CLI adapter schema.
- The full legacy FAVOR annotation schema.
- A new prioritization model.
- A replacement for SeqArray GDS.
- A replacement for IMPACT-VIS input expectations.

## Final Output Summary

The final output of IMPACT-SNV Step 4 is one GDS file per sample:

```text
{sample_id}_SNV_IMPACT.gds
```

Each final GDS file contains the subset of variants retained for that sample after IMPACT prioritization.

At minimum, each final GDS file should contain:

- Core SeqArray genotype and variant nodes.
- Prioritized variants for the sample.
- IMPACT score annotations.
- IMPACT score calculation annotations.
- Tier annotations.
- ClinVar boolean flag annotations.

The output contract must remain stable during the FAVOR-CLI adapter refactor unless an explicit compatibility-breaking change is approved.

## Final File Naming Contract

Final per-sample files must use the naming pattern:

```text
{sample_id}_SNV_IMPACT.gds
```

Examples:

```text
sample_001_SNV_IMPACT.gds
participantA_SNV_IMPACT.gds
```

Do not rename final files to:

```text
{sample_id}.gds
{sample_id}_impact.gds
{sample_id}_SNV.gds
{sample_id}_FAVOR_IMPACT.gds
```

unless a future migration explicitly updates downstream consumers and documentation.

## Final Output Directory Contract

In the current DNAnexus/RAP applet implementation, final output files are collected under the Step 4 output structure and uploaded as the `snv_impact_gds` applet output.

The final applet outputs include:

```text
snv_impact_gds: array of *_SNV_IMPACT.gds files
log_file: impact_prioritization.log
```

Future local, HPC, or container entrypoints may choose different physical output directories, but should preserve the final file naming pattern and the logical output contract.

## Core SeqArray Node Expectations

Final GDS files should remain valid SeqArray-compatible GDS files.

Core nodes expected by downstream tooling may include:

```text
variant.id
chromosome
position
sample.id
genotype
$ref
$alt
```

Not all downstream code may explicitly read every node, but these are part of the expected SeqArray structure and should not be removed by the FAVOR-CLI refactor.

## Sample Contract

Each final GDS file is sample-specific.

Expected behavior:

- The final file name includes the sample identifier.
- The final GDS should contain the relevant sample ID in `sample.id`.
- Sample identity must be preserved from the input GDS.
- Sample-specific variant filtering must not mix variants or genotypes across samples.

The FAVOR-CLI adapter refactor must not change sample identity handling.

## Genotype Preservation Contract

Final GDS files should preserve genotype information for retained variants.

Step 4 currently filters variants per sample based on non-missing genotype information and non-zero IMPACT score.

The refactor must preserve the following principles:

- Genotypes are not recomputed by FAVOR-CLI.
- FAVOR-CLI annotation must not alter genotype calls.
- Adapter injection must not alter genotype nodes.
- Final per-sample files should contain genotype data for retained variants.
- Variants should not be retained or dropped due to annotation row-order mismatch.

Any future change that modifies genotype filtering must be treated as a scientific behavior change.

## Variant Retention Contract

The final `{sample_id}_SNV_IMPACT.gds` files should contain variants that satisfy the current Step 4 retention behavior.

Conceptually, variants are retained for a sample when:

1. The sample has non-missing genotype data for the variant.
2. The variant receives a non-zero IMPACT prioritization score.

The exact filtering implementation may be improved for robustness, but the intended behavior should remain stable unless explicitly changed.

## Required Final IMPACT Annotation Nodes

Each final GDS file should contain the following final IMPACT annotation nodes:

```text
annotation/info/impact_score
annotation/info/impact_score_calc
annotation/info/tier
annotation/info/clnsig_flags/*
```

These nodes are part of the downstream output contract and should not be renamed or relocated during the FAVOR-CLI adapter refactor.

## `impact_score`

### Required path

```text
annotation/info/impact_score
```

### Purpose

Stores the final numeric IMPACT prioritization score for each retained variant.

### Expected type

Numeric vector.

### Expected value range

Conceptually:

```text
0-100
```

Final per-sample GDS files are expected to contain prioritized variants, so retained variants should generally have non-zero scores.

### Source

Current Step 4 initially writes:

```text
annotation/info/patho_score
```

Then post-processing renames or copies this value to:

```text
annotation/info/impact_score
```

### Compatibility rule

Do not rename this node.

Do not replace it with:

```text
patho_score
score
priority_score
favor_score
```

as the final public output field.

## `impact_score_calc`

### Required path

```text
annotation/info/impact_score_calc
```

### Purpose

Stores a character description of the tier and formula used to calculate the IMPACT score.

Examples may include strings conceptually similar to:

```text
Tier 1, 80 + 20 * globalScore
Tier 2, 60 + 40 * globalScore
Tier 3 = 20 + 80 * globalScore
Tier 4 = 100 * ((0.5 * normalized_APC + 0.5 * globalScore))
```

Exact strings may vary based on current Step 4 formatting, but the field must remain interpretable.

### Expected type

Character-like vector.

### Source

Current Step 4 initially writes:

```text
annotation/info/patho_score_calc
```

Then post-processing renames or copies this value to:

```text
annotation/info/impact_score_calc
```

### Compatibility rule

Do not remove this field.

It is useful for transparency, debugging, downstream visualization, and auditability.

## `tier`

### Required path

```text
annotation/info/tier
```

### Purpose

Stores the final IMPACT tier assignment for each retained variant.

### Expected type

Integer vector.

Expected values:

```text
1
2
3
4
```

Missing or unassigned values should be avoided in final retained variant outputs unless explicitly documented.

### Source

Current post-processing extracts tier information from `impact_score_calc` / `patho_score_calc`-style calculation strings.

### Compatibility rule

Do not change tier numbering without an explicit prioritization redesign.

Do not introduce new tiers as part of the initial FAVOR-CLI adapter refactor.

## `clnsig_flags`

### Required path

```text
annotation/info/clnsig_flags/
```

### Purpose

Stores boolean indicator nodes derived from ClinVar clinical significance labels.

These flags support downstream filtering and visualization.

### Required or expected flag nodes

Expected flag nodes include:

```text
annotation/info/clnsig_flags/pathogenic
annotation/info/clnsig_flags/likely_pathogenic
annotation/info/clnsig_flags/uncertain_significance
annotation/info/clnsig_flags/likely_benign
annotation/info/clnsig_flags/benign
annotation/info/clnsig_flags/pathogenic_low_penetrance
annotation/info/clnsig_flags/likely_pathogenic_low_penetrance
annotation/info/clnsig_flags/established_risk_allele
annotation/info/clnsig_flags/likely_risk_allele
annotation/info/clnsig_flags/uncertain_risk_allele
annotation/info/clnsig_flags/affects
annotation/info/clnsig_flags/association
annotation/info/clnsig_flags/drug_response
annotation/info/clnsig_flags/confers_sensitivity
annotation/info/clnsig_flags/protective
annotation/info/clnsig_flags/conflicting_interpretations_of_pathogenicity
annotation/info/clnsig_flags/other
annotation/info/clnsig_flags/not_provided
```

### Expected type

Boolean or bit-packed logical vectors.

### Compatibility rule

Do not remove or rename these flags without updating downstream consumers.

If ClinVar labels are missing from the annotation source, `not_provided` should be populated consistently.

## Prioritization Formula Contract

The final output GDS contract depends on preserving the current IMPACT scoring formulas.

Do not change these formulas during the FAVOR-CLI adapter refactor:

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

Changes to formula definitions are scientific changes and must be handled outside the initial annotation-backend refactor.

## Relationship to Step 4 Intermediate Fields

Current Step 4 may create intermediate fields before final post-processing:

```text
annotation/info/patho_score
annotation/info/patho_score_calc
```

The final public output contract uses:

```text
annotation/info/impact_score
annotation/info/impact_score_calc
```

Implementation details may evolve, but the final output must contain the `impact_*` node names.

## Relationship to Input Annotation Fields

Final output fields depend on the pre-Step-4 annotation compatibility contract.

Before Step 4 runs, the annotated GDS should contain:

```text
annotation/info/FunctionalAnnotation/clnsig
annotation/info/FunctionalAnnotation/genecode_comprehensive_info
annotation/info/FunctionalAnnotation/genecode_comprehensive_exonic_category
annotation/info/FunctionalAnnotation/refseq_exonic_category
annotation/info/FunctionalAnnotation/ucsc_exonic_category
annotation/info/FunctionalAnnotation/apc_protein_function_v3
```

These fields are not necessarily part of the final public output contract, but they are required for producing the final IMPACT annotations.

A future FAVOR-CLI adapter must satisfy the input annotation contract before final output compatibility can be trusted.

## IMPACT-VIS Compatibility

The final `*_SNV_IMPACT.gds` files are intended to be compatible with IMPACT-VIS.

The following should remain stable for IMPACT-VIS compatibility:

- Final file naming pattern: `{sample_id}_SNV_IMPACT.gds`.
- `impact_score` node path.
- `impact_score_calc` node path.
- `tier` node path.
- `clnsig_flags` node paths.
- Core GDS variant/sample/genotype structure.

Any change to these fields may require coordinated updates to IMPACT-VIS documentation and code.

## X/Y Chromosome Output Expectations

The refactor should eventually support final per-sample outputs that include prioritized X and Y chromosome variants when present.

Expected future-compatible intermediate file patterns include:

```text
{sample_id}_chrX.gds
{sample_id}_chrY.gds
```

Expected final merged output remains:

```text
{sample_id}_SNV_IMPACT.gds
```

The final output name should not encode chromosome names because it represents a merged per-sample IMPACT result.

## Provenance Expectations

Final or intermediate GDS files may contain provenance metadata from the FAVOR-CLI adapter and annotation injection steps.

Recommended provenance namespaces:

```text
annotation/info/IMPACT_AnnotationProvenance/
annotation/info/IMPACT_AnnotationCompatibility/
```

Recommended provenance fields include:

```text
favor_cli_version
favor_database_version
reference_genome_build
annotation_command
annotation_date_utc
adapter_version
adapter_schema_version
variant_key_policy
field_mapping_policy
fallback_policy
```

Provenance metadata is strongly recommended for adapter-injected GDS files, especially because updated FAVOR databases may legitimately change annotations and scores.

## Final GDS Validation Requirements

Before a final GDS file is considered valid, validation should confirm:

```text
file exists
file can be opened by SeqArray
sample.id exists
variant.id exists
chromosome exists
position exists
genotype exists
annotation/info/impact_score exists
annotation/info/impact_score_calc exists
annotation/info/tier exists
annotation/info/clnsig_flags/ exists
all final annotation node lengths match variant count
sample identity matches file naming expectation
```

For X/Y support, validation should additionally confirm:

```text
X variants are retained when expected
Y variants are retained when expected
sample-specific chromosome files merge correctly
```

## Score Validation Requirements

Golden validation should include expected score checks for each tier.

Example expectations:

```text
Tier 1 score = 80 + 20 * globalScore
Tier 2 score = 60 + 40 * globalScore
Tier 3 score = 20 + 80 * globalScore
Tier 4 score = 100 * ((0.5 * normalized_APC) + (0.5 * globalScore))
```

Validation should distinguish:

- Formula regressions.
- Annotation-content changes due to updated FAVOR data.
- Adapter mapping errors.
- Variant filtering differences.

## Backward Compatibility Requirements

During the FAVOR-CLI refactor:

- Existing legacy FAVORannotator-produced annotated GDS files should still be accepted by Step 4.
- Existing final output naming must remain unchanged.
- Existing final IMPACT node names must remain unchanged.
- Existing `GeneList.txt` format must remain unchanged.
- Existing IMPACT-VIS expectations must remain valid.

## Compatibility-Breaking Changes

The following are compatibility-breaking and should not be done without explicit approval:

- Renaming `{sample_id}_SNV_IMPACT.gds` outputs.
- Removing `annotation/info/impact_score`.
- Removing `annotation/info/impact_score_calc`.
- Removing `annotation/info/tier`.
- Removing or renaming `annotation/info/clnsig_flags/*` nodes.
- Changing Step 4 scoring formulas.
- Changing `GeneList.txt` required columns.
- Replacing GDS as the final output format.
- Changing sample-specific output semantics.

## Testing Expectations

Output GDS tests should include:

- A final GDS with Tier 1 variant.
- A final GDS with Tier 2 variant.
- A final GDS with Tier 3 variant.
- A final GDS with Tier 4 variant.
- A final GDS with multiple samples processed separately.
- A final GDS containing an X chromosome variant.
- A final GDS containing a Y chromosome variant.
- A final GDS with missing ClinVar source annotation but valid fallback flags.
- A final GDS with expected `not_provided` ClinVar flag behavior.

## Implementation Guidance for Coding Agents

When modifying code that affects final outputs:

- Do not change final output names.
- Do not change final IMPACT node paths.
- Do not change Step 4 formulas.
- Do not change `GeneList.txt` format.
- Add validation before changing output generation.
- Confirm final GDS files can be opened with SeqArray.
- Confirm final annotations have lengths equal to variant count.
- Summarize any output-impacting behavior change explicitly.

## Future Revisions

This contract may be revised only when the output format is intentionally changed.

Any future revision should document:

- Why the output contract changed.
- Whether IMPACT-VIS was updated.
- Whether migration tools are needed.
- Whether old outputs remain readable.
- Whether final file naming changed.
- Whether score interpretation changed.

## Final Principle

The FAVOR-CLI refactor may change annotation sources and annotation freshness.

It should not silently change the final IMPACT-SNV output contract.

The final `{sample_id}_SNV_IMPACT.gds` files must remain structurally and semantically compatible with current downstream expectations unless a deliberate compatibility-breaking migration is approved.
