# Copilot Instructions for IMPACT-SNV

IMPACT-SNV is a computational genomics pipeline for SNV/indel processing, FAVOR annotation, variant prioritization, and generation of per-sample GDS files for downstream IMPACT-VIS visualization.

## Primary Refactor Goal

We are refactoring IMPACT-SNV to replace the legacy FAVORannotator-based annotation step with a future FAVOR-CLI-backed annotation adapter.

The intended architecture is:

```text
VCF/GDS variant identity
→ FAVOR-CLI annotation
→ normalized IMPACT annotation compatibility schema
→ inject annotations into SeqArray GDS
→ existing IMPACT prioritization
→ final {sample_id}_SNV_IMPACT.gds
```

The refactor should modernize annotation while preserving the existing IMPACT-SNV scientific and output contracts.

## Current Pipeline Context

The current repository is organized as a four-step workflow:

```text
VCF files
→ step1_vcf_merge/
→ step2_vcf2gds/
→ Step3_favorannotator-rap/
→ step4_impact_prioritization/
→ final {sample_id}_SNV_IMPACT.gds
```

Current responsibilities:

- `step1_vcf_merge/`: merges input VCFs, keeps GT FORMAT, normalizes variants, and splits by chromosome.
- `step2_vcf2gds/`: converts chromosome-specific VCF files to SeqArray GDS and adds `annotation/info/QC_label = PASS`.
- `Step3_favorannotator-rap/`: legacy FAVORannotator-based annotation step that writes FAVOR-style annotations into GDS.
- `step4_impact_prioritization/`: reads annotated GDS files, applies IMPACT prioritization, exports per-sample GDS files, and post-processes final IMPACT annotations.

## Critical Guardrails

Do not change the following unless explicitly requested:

- Step 4 scoring formulas.
- `GeneList.txt` input contract.
- Final output naming: `{sample_id}_SNV_IMPACT.gds`.
- Final output annotation nodes:
  - `annotation/info/impact_score`
  - `annotation/info/impact_score_calc`
  - `annotation/info/tier`
  - `annotation/info/clnsig_flags/*`
- SeqArray GDS as the internal downstream exchange format.
- Legacy `Step3_favorannotator-rap/` behavior.
- Compatibility with existing legacy FAVORannotator-produced annotated GDS files.

Do not delete or rewrite legacy FAVORannotator code until the FAVOR-CLI adapter path has been validated.

## Current Step 4 Required Annotation Nodes

Step 4 currently expects legacy-style FAVOR annotations under:

```text
annotation/info/FunctionalAnnotation/
```

Required nodes include:

```text
annotation/info/FunctionalAnnotation/clnsig
annotation/info/FunctionalAnnotation/genecode_comprehensive_info
annotation/info/FunctionalAnnotation/genecode_comprehensive_exonic_category
annotation/info/FunctionalAnnotation/refseq_exonic_category
annotation/info/FunctionalAnnotation/ucsc_exonic_category
annotation/info/FunctionalAnnotation/apc_protein_function_v3
```

Any FAVOR-CLI adapter must provide these fields, or a clearly validated compatibility equivalent, before Step 4 runs.

## Prioritization Formulas Must Remain Stable

Do not change the IMPACT prioritization formulas unless explicitly requested.

Current formulas:

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

`GeneList.txt` must continue to use:

```text
symbol
globalScore
```

## FAVOR-CLI Refactor Principles

When implementing FAVOR-CLI support:

- Treat FAVOR-CLI as an annotation backend, not as a replacement for the full IMPACT-SNV pipeline.
- Prefer adding a new `step3_favorcli_annotation/` module rather than modifying legacy `Step3_favorannotator-rap/` in place.
- Do not assume the exact FAVOR-CLI output schema unless real output examples or official schema output are available.
- If schema details are uncertain, create explicit placeholders, fixtures, mapping configuration, or validation hooks.
- Normalize FAVOR-CLI output into an IMPACT-owned compatibility schema before writing to GDS.
- Fail fast when required fields are missing or unsafe to map.
- Record all fallbacks explicitly.
- Preserve source/provenance metadata whenever annotations are generated or injected.

## GDS Annotation Injection Requirements

Adapter-injected annotations should be written into the same namespace consumed by current Step 4:

```text
annotation/info/FunctionalAnnotation/
```

Required compatibility fields should be readable with `seqGetData()`.

Before Step 4 runs, validate:

- Required nodes exist.
- Required nodes have expected types.
- Required node lengths match the GDS variant count.
- Variant keys and row ordering are aligned.
- Fallbacks and missing fields are reported.

Store provenance under stable namespaces such as:

```text
annotation/info/IMPACT_AnnotationProvenance/
annotation/info/IMPACT_AnnotationCompatibility/
```

Useful provenance fields include:

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

## Chromosome Handling

The refactor should support chromosomes:

```text
1-22, X, Y
```

Preserve filename conventions such as:

```text
merged_chr1.vcf.gz
merged_chrX.vcf.gz
merged_chrY.vcf.gz
merged_chr1.gds
merged_chrX.gds
merged_chrY.gds
{sample_id}_chrX.gds
{sample_id}_chrY.gds
```

Do not introduce mitochondrial support as required behavior unless explicitly requested. If added, make it optional and clearly documented.

Normalize chromosome aliases carefully:

```text
chr1 <-> 1
chrX <-> X
chrY <-> Y
```

Do not silently convert `23` to `X` or `24` to `Y` unless explicitly configured.

## Implementation Style

Use small, reviewable PR-sized changes.

Preferred implementation order:

```text
1. Documentation and contracts
2. Chromosome handling generalization for 1-22, X, Y
3. FAVOR-CLI Step 3 skeleton with dry-run mode
4. Annotation schema adapter using real or fixture FAVOR-CLI output
5. GDS annotation injection
6. Minimal Step 4 robustness updates
7. Golden-dataset validation
8. Local/HPC/container/RAP portability improvements
9. Documentation and migration guide
10. Legacy FAVORannotator deprecation gate
```

For each task:

- Inspect relevant files before editing.
- Preserve backward compatibility.
- Avoid broad rewrites.
- Add validation before changing behavior.
- Prefer explicit configuration over hard-coded assumptions.
- Keep RAP/DNAnexus support, but separate platform wrappers from core logic where practical.
- Do not introduce a full workflow engine unless explicitly requested.

## Testing and Validation Expectations

When adding behavior, include at least one of:

- Unit tests.
- Small fixture-based tests.
- Dry-run validation.
- Preflight schema validation.
- Documented manual validation commands.

Golden validation should eventually check:

- Required GDS nodes exist.
- Node lengths match variant count.
- Sample IDs are preserved.
- Genotypes are preserved for retained variants.
- X/Y variants are not dropped.
- Scores and tiers match expected formulas.
- Final `{sample_id}_SNV_IMPACT.gds` outputs are produced.

## Do Not Refactor Yet

Do not implement these changes unless explicitly requested:

- Do not rewrite Step 4 scoring logic.
- Do not replace GDS with Parquet or VCF as the downstream format.
- Do not introduce Nextflow, Snakemake, WDL, or another workflow engine as a required dependency.
- Do not redesign the annotation model around new FAVOR features, tissue-specific annotations, ACMG-like rules, or new prioritization tiers.
- Do not remove legacy FAVORannotator support.
- Do not silently synthesize scientifically meaningful annotations without reporting fallbacks.

## Communication Expectations

When completing a task, summarize:

- Files changed.
- Behavior changed, if any.
- Backward compatibility impact.
- Tests or validation performed.
- Known limitations or follow-up tasks.

If a requested change could alter scientific interpretation, call that out explicitly before implementing it.
