# IMPACT-SNV Documentation

This directory contains comprehensive documentation for the IMPACT-SNV pipeline and its ongoing refactor to modernize annotation infrastructure.

## Overview

The IMPACT-SNV pipeline is organized as four sequential processing steps that transform raw VCF files into prioritized per-sample GDS files compatible with [IMPACT-VIS](https://boehlernick.github.io/IMPACT-VIS).

All documentation in this directory is designed to support developers, contributors, and automated coding agents during implementation and refactoring work.

## Documentation Files

### Pipeline Reference

#### [pipeline_contract.md](pipeline_contract.md)
**Current state specification of the IMPACT-SNV pipeline**

Describes the existing four-step workflow before refactoring:
- Step 1: VCF merge, normalization, and chromosome split
- Step 2: VCF-to-GDS conversion
- Step 3: Legacy FAVORannotator annotation
- Step 4: IMPACT variant prioritization

Specifies inputs, outputs, file naming conventions, and known behavior for each step. This is the authoritative reference for current pipeline behavior.

### Refactoring Guidance

#### [favorcli_refactor_plan.md](favorcli_refactor_plan.md)
**Staged refactor strategy to modernize annotation infrastructure**

Outlines the high-level refactoring approach:
- Separation of FAVOR-CLI as a replaceable annotation backend
- Preservation of IMPACT prioritization logic and scoring methods
- Phased migration strategy
- Architecture for adapter layer between FAVOR-CLI and IMPACT-SNV

Read this first to understand the refactoring vision and rationale.

### Contract Specifications

#### [annotation_compatibility_contract.md](annotation_compatibility_contract.md)
**Required annotation fields and compatibility schema for IMPACT prioritization**

Defines the minimum annotation interface required by Step 4:
- Required GDS annotation node paths and naming
- Expected annotation types and value conventions
- Current annotation consumers in Step 4 prioritization
- Normalization requirements for new annotation sources

**Critical reference** for any new annotation adapter (FAVOR-CLI or otherwise) that must supply variants to Step 4.

#### [output_gds_contract.md](output_gds_contract.md)
**Final output file structure and preservation requirements**

Specifies the structure and stability guarantees for final `{sample_id}_SNV_IMPACT.gds` files:
- File naming pattern
- Core SeqArray node expectations
- IMPACT score and tier annotations
- ClinVar boolean flag annotations
- Sample and genotype preservation requirements

Ensures downstream tools (IMPACT-VIS, etc.) receive compatible outputs during refactoring.

#### [chromosome_handling_contract.md](chromosome_handling_contract.md)
**Current chromosome handling behavior and refactor targets**

Documents chromosome support across the pipeline:
- Current limitation: autosome-only (chr1-chr22)
- Target support: chr1-chr22, chrX, chrY
- Chromosome naming conventions
- Validation expectations for X/Y expansion

Provides context for planned chromosome support improvements.

### Discovery and Planning

#### [favorcli_schema_discovery.md](favorcli_schema_discovery.md)
**Protocol for empirically discovering FAVOR-CLI output schema**

Guidelines for schema discovery before implementing the FAVOR-CLI adapter:
- How to run FAVOR-CLI on controlled test data
- How to inspect and document output schema
- How to compare FAVOR-CLI output against required annotation compatibility schema
- Documentation requirements for discovered mappings

Prevents guessing about FAVOR-CLI field mappings and ensures adapter decisions are evidence-based.

## Getting Started

### For Implementation Work
1. Start with **[pipeline_contract.md](pipeline_contract.md)** to understand current behavior
2. Read **[favorcli_refactor_plan.md](favorcli_refactor_plan.md)** for refactoring strategy
3. Review the relevant **contract** document for your work:
   - Working on annotation: **[annotation_compatibility_contract.md](annotation_compatibility_contract.md)**
   - Working on output: **[output_gds_contract.md](output_gds_contract.md)**
   - Working on chromosomes: **[chromosome_handling_contract.md](chromosome_handling_contract.md)**

### For FAVOR-CLI Adapter Development
1. Complete schema discovery using **[favorcli_schema_discovery.md](favorcli_schema_discovery.md)**
2. Review **[annotation_compatibility_contract.md](annotation_compatibility_contract.md)** for required fields
3. Document field mappings and transformations
4. Implement adapter with reference to **[pipeline_contract.md](pipeline_contract.md)** Step 3

### For Code Review
1. Reference **[pipeline_contract.md](pipeline_contract.md)** for expected behavior
2. Check **[annotation_compatibility_contract.md](annotation_compatibility_contract.md)** for annotation requirements
3. Verify **[output_gds_contract.md](output_gds_contract.md)** for output stability
4. Consult the relevant **contract** document if behavior changes are proposed

## Key Concepts

### Annotation Compatibility Schema
The minimum set of GDS annotation nodes required by IMPACT-SNV Step 4 prioritization:
```
annotation/info/FunctionalAnnotation/clnsig
annotation/info/FunctionalAnnotation/genecode_comprehensive_info
annotation/info/FunctionalAnnotation/genecode_comprehensive_exonic_category
annotation/info/FunctionalAnnotation/refseq_exonic_category
annotation/info/FunctionalAnnotation/ucsc_exonic_category
annotation/info/FunctionalAnnotation/apc_protein_function_v3
```

Any new annotation source (including FAVOR-CLI) must provide these fields or validated equivalents.

### Prioritization Formulas
Step 4 assigns pathogenicity scores using tier-based formulas:
- **Tier 1**: ClinVar Pathogenic/Likely Pathogenic + gene match = `80 + 20 × globalScore`
- **Tier 2**: Frameshift/stopgain = `60 + 40 × globalScore`
- **Tier 3**: Nonsynonymous/nonframeshift = `20 + 80 × globalScore`
- **Tier 4**: APC protein function = `100 × (0.5 × normalized_APC + 0.5 × globalScore)`

These formulas must remain stable during refactoring unless explicitly changed.

### Output Format
Final output files follow the pattern: `{sample_id}_SNV_IMPACT.gds`

Each file contains one sample's prioritized variants with IMPACT scores, tiers, and ClinVar flags.

## Additional Resources

- **Root Repository README**: [../README.md](../README.md) - Pipeline overview and usage
- **Individual Step READMEs**:
  - [../step1_vcf_merge/README.md](../step1_vcf_merge/README.md) or README.developer.md
  - [../step2_vcf2gds/README.md](../step2_vcf2gds/README.md) or README.app.md
  - [../Step3_favorannotator-rap/README.md](../Step3_favorannotator-rap/README.md) or Readme.app.md
  - [../step4_impact_prioritization/README.md](../step4_impact_prioritization/README.md)
- **GitHub Copilot Instructions**: [../.github/copilot-instructions.md](../.github/copilot-instructions.md)
- **IMPACT-VIS**: [https://github.com/boehlernick/IMPACT-VIS](https://github.com/boehlernick/IMPACT-VIS)

## Maintenance

These documents should be updated when:
- Pipeline behavior changes
- New annotation fields are introduced
- Output structure is modified
- New chromosome support is added
- FAVOR-CLI adapter is implemented

All changes to pipeline behavior, annotation requirements, or output structure should include documentation updates to maintain this reference accuracy.
