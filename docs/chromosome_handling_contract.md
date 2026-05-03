# IMPACT-SNV Chromosome Handling Contract

## Purpose

This document defines the chromosome handling contract for IMPACT-SNV before and during the FAVOR-CLI refactor.

It exists to support a safe expansion from the current autosome-oriented workflow to a workflow that supports:

```text
1-22, X, Y
```

The contract defines chromosome naming, normalization, expected file naming conventions, validation expectations, and behavior that must remain stable during refactoring.

This document should be read alongside:

- `.github/copilot-instructions.md`
- `docs/favorcli_refactor_plan.md`
- `docs/pipeline_contract.md`
- `docs/annotation_compatibility_contract.md`
- `docs/output_gds_contract.md`
- `docs/favorcli_schema_discovery.md`

## Scope

This document covers:

- Current chromosome handling behavior.
- Target chromosome support for the refactor.
- Chromosome naming conventions.
- Chromosome alias normalization.
- Step-specific expectations for Steps 1-4.
- X and Y chromosome handling.
- Optional future mitochondrial support.
- Validation and testing expectations.

This document does **not** define:

- Sex-aware genotype interpretation.
- PAR-region-specific biological interpretation.
- Ploidy-aware scoring changes.
- New prioritization formulas for sex chromosomes.
- Required mitochondrial chromosome support.

The initial refactor goal is to stop dropping X and Y variants and to allow them to flow through the same IMPACT-SNV prioritization framework where annotation support exists.

The current repository state already includes X/Y-aware chromosome splitting in Step 1, so the remaining work in this contract is mainly to preserve that behavior and avoid accidental regressions.

## Current Behavior Summary

The current IMPACT-SNV workflow is primarily autosome-oriented.

Current Step 1 behavior:

- Splits merged VCF output into chromosome-specific files.
- Attempts chromosome-specific outputs for `chr1` through `chr22`, `chrX`, and `chrY`.
- Preserves the chromosome naming style already present in the merged VCF when choosing `bcftools view --regions` inputs.
- Skips chromosomes with no variants instead of failing the step.

Current Step 2 behavior:

- Converts provided VCF files to GDS.
- Does not inherently need to be autosome-only if X/Y VCF files are provided.

Current Step 3 legacy FAVORannotator behavior:

- Historically expects autosome-style chromosome input.
- Uses legacy FAVOR database assumptions.
- Should not be relied upon for X/Y expansion during the FAVOR-CLI refactor.

Current Step 4 behavior:

- Already contains filename patterns that can match `X`, `Y`, or `M` chromosome suffixes.
- Can discover `merged_chrX.gds` and `merged_chrY.gds` when present.
- Still depends on earlier steps producing and annotating those chromosome files.

## Target Chromosome Support

The initial target chromosome set for the FAVOR-CLI refactor is:

```text
1, 2, 3, ..., 22, X, Y
```

Equivalent display/output form:

```text
chr1, chr2, chr3, ..., chr22, chrX, chrY
```

Mitochondrial support is not required in the initial refactor.

If mitochondrial support is added later, it must be optional and explicitly documented.

## Canonical Internal Chromosome Representation

Recommended canonical internal representation:

```text
1
2
...
22
X
Y
```

That is, chromosome values should be normalized internally without the `chr` prefix.

File naming may still use the `chr` prefix for backward compatibility.

## Accepted External Aliases

The pipeline should accept common chromosome aliases where practical.

Recommended accepted aliases:

```text
chr1  -> 1
chr2  -> 2
...
chr22 -> 22
chrX  -> X
chrY  -> Y
1     -> 1
2     -> 2
...
22    -> 22
X     -> X
Y     -> Y
```

Case normalization should be supported for sex chromosomes:

```text
x -> X
y -> Y
chrx -> X
chry -> Y
```

However, output file names should preserve the expected `chr`-prefixed convention.

## Aliases That Must Not Be Silently Converted

Do not silently convert numeric sex chromosome aliases such as:

```text
23 -> X
24 -> Y
```

unless that behavior is explicitly configured, documented, and recorded in provenance.

Rationale:

- Different datasets and tools may encode sex chromosomes differently.
- Silent conversion can hide input inconsistencies.
- Scientific interpretation and file matching can become ambiguous.

If a future implementation supports these aliases, it should require an explicit option such as:

```text
allow_numeric_sex_chromosomes: true
```

and should record the normalization policy in logs or provenance metadata.

## Unsupported or Optional Chromosomes

The initial refactor does not require support for:

```text
M
MT
chrM
chrMT
mitochondrial
```

If mitochondrial support is added later, accepted aliases may include:

```text
M
MT
chrM
chrMT
```

but mitochondrial support must remain optional unless future documentation explicitly updates this contract.

Unplaced contigs, alternate haplotypes, decoy contigs, and patch sequences are out of scope for the initial FAVOR-CLI refactor unless explicitly requested.

Examples out of scope:

```text
GL000207.1
KI270728.1
chrUn
chr1_KI270706v1_random
HLA contigs
alternate haplotypes
```

## File Naming Contract

The current and future-compatible pipeline should use the following chromosome-specific file naming patterns.

### Step 1 VCF outputs

Autosomes:

```text
merged_chr1.vcf.gz
merged_chr2.vcf.gz
...
merged_chr22.vcf.gz
```

Sex chromosomes:

```text
merged_chrX.vcf.gz
merged_chrY.vcf.gz
```

### Step 2 GDS outputs

Autosomes:

```text
merged_chr1.gds
merged_chr2.gds
...
merged_chr22.gds
```

Sex chromosomes:

```text
merged_chrX.gds
merged_chrY.gds
```

### Step 4 sample-specific chromosome intermediate outputs

Autosomes:

```text
{sample_id}_chr1.gds
{sample_id}_chr2.gds
...
{sample_id}_chr22.gds
```

Sex chromosomes:

```text
{sample_id}_chrX.gds
{sample_id}_chrY.gds
```

### Final merged per-sample output

The final output must remain:

```text
{sample_id}_SNV_IMPACT.gds
```

Do not include chromosome names in the final merged per-sample output file name.

## Step 1 Contract: VCF Splitting

Step 1 should eventually support splitting merged VCF files into:

```text
merged_chr1.vcf.gz
...
merged_chr22.vcf.gz
merged_chrX.vcf.gz
merged_chrY.vcf.gz
```

### Required behavior

Step 1 should:

- Preserve existing autosome behavior.
- Include X and Y when present in input data.
- Avoid dropping X/Y variants silently.
- Use consistent output naming.
- Log which chromosomes were requested.
- Log which chromosome outputs were produced.
- Log which requested chromosomes had zero variants or no output.

### Empty chromosome behavior

The implementation must explicitly choose and document one behavior:

1. Emit empty chromosome VCF files, or
2. Skip chromosomes with no variants.

Recommended initial behavior:

```text
Skip chromosome-specific outputs for chromosomes with no variants, but log the skipped chromosomes clearly.
```

Rationale:

- Empty VCF/GDS files can cause unnecessary downstream failures.
- Skipping empty chromosomes is often easier to handle as long as it is explicit.

If empty files are emitted, downstream steps must be able to handle them safely.

### Input chromosome naming

Step 1 must account for whether the input VCF uses:

```text
chr1, chr2, ..., chrX, chrY
```

or:

```text
1, 2, ..., X, Y
```

The splitting logic should not assume `chr` prefixes unless the input has been normalized to that convention.

Recommended behavior:

- Detect chromosome naming style from the VCF header or records.
- Use matching region names when calling `bcftools view --regions`.
- Preserve output file naming as `merged_chr*.vcf.gz` regardless of internal input style.

## Step 2 Contract: VCF-to-GDS Conversion

Step 2 should treat X/Y VCF files the same way it treats autosomal VCF files where possible.

Expected input examples:

```text
merged_chr1.vcf.gz
merged_chrX.vcf.gz
merged_chrY.vcf.gz
```

Expected output examples:

```text
merged_chr1.gds
merged_chrX.gds
merged_chrY.gds
```

Step 2 should not introduce chromosome-specific scoring logic.

If SeqArray handles X/Y chromosome labels differently based on input naming, this behavior should be documented and validated.

## Step 3 Contract: Annotation

### Legacy Step 3

The legacy FAVORannotator-based Step 3 should not be modified to claim full X/Y support unless explicitly implemented and validated.

The initial X/Y-supporting annotation path should be developed through the future FAVOR-CLI-backed annotation adapter.

### Future FAVOR-CLI Step 3

The FAVOR-CLI annotation path should support:

```text
1-22, X, Y
```

The adapter must ensure that chromosome labels in FAVOR-CLI output align with chromosome labels in the GDS variant identity table.

The adapter must validate:

```text
chromosome labels match
positions match
REF alleles match
ALT alleles match
variant row order is preserved or correctly restored
```

If FAVOR-CLI uses one naming convention and GDS uses another, normalization must be explicit.

Example:

```text
FAVOR-CLI output: chrX
GDS canonical key: X-12345-A-G
normalized key: X-12345-A-G
```

## Step 4 Contract: Prioritization and Sample Merge

Step 4 should process X/Y chromosome GDS files when they are present and annotated.

Expected input examples:

```text
merged_chr1.gds
merged_chr2.gds
merged_chrX.gds
merged_chrY.gds
```

Expected sample-specific intermediate outputs:

```text
{sample_id}_chr1.gds
{sample_id}_chrX.gds
{sample_id}_chrY.gds
```

Expected final output:

```text
{sample_id}_SNV_IMPACT.gds
```

### Required behavior

Step 4 should:

- Discover X/Y chromosome GDS files when present.
- Extract chromosome labels correctly from filenames.
- Export sample-specific X/Y intermediate files where variants are retained.
- Merge sample-specific X/Y files into final per-sample GDS outputs.
- Preserve existing autosome behavior.
- Preserve current scoring formulas.

### Forbidden behavior

Step 4 must not:

- Drop X/Y files because of autosome-only regex patterns.
- Change scoring formulas for X/Y variants.
- Rename final output files to include chromosome names.
- Require all samples to have Y chromosome variants.
- Fail solely because a chromosome has no retained variants, if the absence is expected.

## Ploidy and Sex Chromosome Considerations

The initial X/Y refactor does not introduce sex-aware or ploidy-aware scoring changes.

Known considerations:

- Male samples may have hemizygous genotypes on parts of chrX or chrY.
- Female samples may have no chrY variants.
- PAR regions may have different biological interpretation.
- Some VCFs may encode sex chromosome genotypes differently.

Initial contract:

- Preserve genotype information as represented in the input GDS.
- Do not reinterpret sex chromosome genotypes.
- Do not change scoring formulas based on sex or ploidy.
- Do not require sample sex metadata for the initial X/Y support.

Future sex-aware interpretation may be valuable, but it is outside the initial FAVOR-CLI refactor.

## PAR Region Considerations

Pseudoautosomal regions are biologically important but are not part of the initial chromosome-handling implementation contract.

Initial behavior:

- PAR variants should be retained if present, annotated, scored, and genotype criteria are met.
- No PAR-specific scoring changes should be introduced.
- No PAR-specific duplicate handling should be added unless required by observed input behavior.

If PAR-specific handling is added later, it must be documented as a scientific behavior change.

## Y Chromosome Missingness

The pipeline must tolerate missing chrY data.

Expected valid scenarios:

- No `merged_chrY.vcf.gz` because no chrY variants exist.
- No `merged_chrY.gds` because no chrY VCF was produced.
- No `{sample_id}_chrY.gds` because a sample has no retained Y variants.
- Final `{sample_id}_SNV_IMPACT.gds` still exists if the sample has retained variants on other chromosomes.

The absence of chrY should not be treated as an error unless chrY was explicitly required for a specific run.

## Chromosome Configuration

Future implementations may define chromosome behavior in a configuration file.

Possible configuration file:

```text
config/chromosomes.yml
```

Possible contents:

```yaml
chromosomes:
  default:
    - "1"
    - "2"
    - "3"
    - "4"
    - "5"
    - "6"
    - "7"
    - "8"
    - "9"
    - "10"
    - "11"
    - "12"
    - "13"
    - "14"
    - "15"
    - "16"
    - "17"
    - "18"
    - "19"
    - "20"
    - "21"
    - "22"
    - "X"
    - "Y"
  include_mitochondrial: false
  allow_numeric_sex_chromosomes: false
  mitochondrial_aliases:
    - "M"
    - "MT"
    - "chrM"
    - "chrMT"
```

Configuration should not be required for the first implementation if a simple constant is sufficient, but chromosome behavior should be explicit and testable.

## Validation Requirements

Chromosome handling changes should validate:

```text
chr1 and 1 normalize consistently
chr22 and 22 normalize consistently
chrX, X, chrx, and x normalize to X
chrY, Y, chry, and y normalize to Y
23 is not silently converted to X
24 is not silently converted to Y
unsupported contigs are reported or skipped explicitly
```

Pipeline-level validation should confirm:

```text
Step 1 can produce merged_chrX.vcf.gz when X variants are present
Step 1 can produce merged_chrY.vcf.gz when Y variants are present
Step 2 can convert merged_chrX.vcf.gz to merged_chrX.gds
Step 2 can convert merged_chrY.vcf.gz to merged_chrY.gds
Step 4 can process merged_chrX.gds
Step 4 can process merged_chrY.gds
Step 4 can merge sample-specific X/Y GDS files into final *_SNV_IMPACT.gds
```

## Testing Expectations

Tests or fixtures should include:

- Autosomal-only input.
- Input containing chrX variants.
- Input containing chrY variants.
- Input containing chr1, chrX, and chrY variants together.
- Input using `chr`-prefixed chromosome names.
- Input using non-prefixed chromosome names.
- A sample with no Y variants.
- A chromosome requested but absent from input.
- An unsupported contig.

Expected test assertions:

```text
X variants are not dropped
Y variants are not dropped
Autosome behavior is unchanged
File naming remains stable
Final per-sample outputs remain *_SNV_IMPACT.gds
No scoring formulas are changed
```

## Logging Requirements

Chromosome-aware steps should log:

```text
requested chromosome set
input chromosome naming style
normalization policy
chromosomes found in input
chromosomes emitted
chromosomes skipped due to no variants
unsupported contigs encountered
```

For FAVOR-CLI adapter runs, chromosome normalization policy should also be recorded in provenance metadata.

## Provenance Requirements

When annotations are generated or injected, record chromosome normalization information where practical.

Recommended provenance fields:

```text
chromosome_normalization_policy
requested_chromosomes
observed_chromosomes
emitted_chromosomes
skipped_chromosomes
unsupported_chromosomes
allow_numeric_sex_chromosomes
include_mitochondrial
```

Suggested GDS namespace:

```text
annotation/info/IMPACT_AnnotationProvenance/
```

## Backward Compatibility Requirements

The chromosome refactor must preserve:

- Existing autosome-only workflows.
- Existing `merged_chr1` through `merged_chr22` naming.
- Existing Step 2 behavior for autosomal VCFs.
- Existing Step 4 scoring formulas.
- Existing final output naming.
- Existing legacy FAVORannotator compatibility for autosomes.

Adding X/Y support should be additive and should not break autosome processing.

## Compatibility-Breaking Changes

The following are compatibility-breaking and should not be done without explicit approval:

- Renaming autosomal files from `merged_chr1` to `merged_1`.
- Removing support for `chr`-prefixed filenames.
- Requiring Y chromosome data for all samples.
- Changing final output naming to include chromosome labels.
- Changing Step 4 scoring based on chromosome.
- Converting `23` to `X` or `24` to `Y` silently.
- Making mitochondrial support mandatory without validation.

## Implementation Guidance for Coding Agents

When implementing chromosome changes:

- Start with Step 1 chromosome splitting.
- Preserve autosome behavior first.
- Add X/Y support additively.
- Avoid broad rewrites.
- Keep file naming stable.
- Do not change Step 4 formulas.
- Add or document validation.
- Make empty chromosome behavior explicit.
- Do not modify legacy FAVORannotator unless explicitly requested.

Recommended first behavior-change task:

```text
Generalize Step 1 chromosome splitting from chr1- chr22 to chr1- chr22 plus chrX and chrY, preserving existing autosome output naming and documenting empty chromosome behavior.
```

Recommended second behavior-change task:

```text
Verify and minimally update Step 4 filename discovery and sample merge logic so merged_chrX.gds, merged_chrY.gds, {sample_id}_chrX.gds, and {sample_id}_chrY.gds are handled when present.
```

## Future Revisions

This contract may be revised if future work adds:

- Mitochondrial support.
- PAR-specific handling.
- Sex-aware genotype interpretation.
- Ploidy-aware scoring.
- Alternate contig support.
- Workflow-engine-specific chromosome configuration.

Any future revision must document whether behavior changes are compatibility-preserving or compatibility-breaking.

## Final Principle

The chromosome handling refactor should be additive and conservative.

The immediate goal is to ensure X and Y variants can flow through the IMPACT-SNV pipeline without being silently dropped, while preserving current autosomal behavior and final output compatibility.
