# IMPACT-SNV

**Integrated Mapping of Phenotype-Associated Candidate Targets for SNV/Indel Analysis**

This repository contains the IMPACT-SNV pipeline, which processes and prioritizes single nucleotide variants (SNVs) and indels for rare disease analysis using the FAVOR database and phenotype-specific gene-disease associations.

## Overview

IMPACT-SNV is part of the broader IMPACT framework for phenotype-configurable interpretation of genomic variants. This module specifically handles SNV/Indel processing and produces output files compatible with [IMPACT-VIS](https://boehlernick.github.io/IMPACT-VIS) for interactive visualization and analysis.

## Pipeline Architecture

The pipeline consists of four modular steps, each available as a DNAnexus applet:

```
VCF Files → [Step 1: Merge] → [Step 2: VCF2GDS] → [Step 3: FAVOR Annotate] → [Step 4: Prioritize] → *_SNV_IMPACT.gds
```

| Step | Folder | Description | Input | Output |
|------|--------|-------------|-------|--------|
| 1 | `step1_vcf_merge/` | Merges multiple VCF files into chromosome-separated files | `.vcf`, `.vcf.gz` | `merged_chr*.vcf.gz` |
| 2 | `step2_vcf2gds/` | Converts VCF to GDS format for efficient processing | `.vcf.gz` | `merged_chr*.gds` |
| 3 | `step3_favorannotator-rap/` | Annotates variants using FAVOR database | `.gds` | `favor_merged_chr*.gds` |
| 4 | `step4_impact_prioritization/` | Scores and prioritizes variants by pathogenicity | `.gds`, `GeneList.txt` | `*_SNV_IMPACT.gds` |

## Requirements

### Software Dependencies
- **R** ≥ 4.0 with Bioconductor packages:
  - `SeqArray`, `SeqVarTools`, `gdsfmt`
  - `stringi`, `dplyr`, `tidyr`
- **Python 3** (for VCF merging scripts)
- **BCFtools**, **bgzip**, **tabix** (for VCF processing)
- **Rust** and **xsv** (for FAVOR annotation step)

### Platform
These tools are designed for deployment on the **DNAnexus Research Analysis Platform**. Each folder contains a `dxapp.json` configuration file for building DNAnexus applets.

## Input Requirements

### Step 1: VCF Merge
- **Input**: One or more VCF or VCF.GZ files
- **Output**: Chromosome-separated VCF.GZ files (`merged_chr1.vcf.gz`, `merged_chr2.vcf.gz`, etc.)

### Step 2: VCF to GDS
- **Input**: VCF.GZ file from Step 1
- **Output**: GDS file (`merged_chr*.gds`)

### Step 3: FAVOR Annotation
- **Input**: GDS file from Step 2
- **Output**: Annotated GDS with FAVOR functional annotations

### Step 4: IMPACT Prioritization
- **Input**:
  - Annotated GDS files from Step 3 (`merged_chr*.gds`)
  - Gene-disease association file (`GeneList.txt`) - tab-separated with columns:
    - `symbol`: Gene symbol (e.g., `GJB2`, `OTOF`)
    - `globalScore`: Open Targets association score (0-1)
- **Output**: Per-sample GDS files (`{sample_id}_SNV_IMPACT.gds`)

## Output File Format

The final `*_SNV_IMPACT.gds` files contain:

### Core SeqArray Nodes
| Node | Type | Description |
|------|------|-------------|
| `variant.id` | integer | Unique variant identifier |
| `chromosome` | integer/character | Chromosome |
| `position` | integer | 1-based genomic position |
| `sample.id` | character | Sample identifiers |
| `genotype` | integer | Genotype array |

### IMPACT Score Annotations
| Node | Type | Description |
|------|------|-------------|
| `annotation/info/impact_score` | numeric | Pathogenicity score (0-100) |
| `annotation/info/impact_score_calc` | character | Tier and calculation formula |
| `annotation/info/tier` | integer | Priority tier (1-4) |

### ClinVar Significance Flags
Boolean indicators under `annotation/info/clnsig_flags/`:
- `pathogenic`, `likely_pathogenic`, `uncertain_significance`
- `likely_benign`, `benign`
- `pathogenic_low_penetrance`, `likely_pathogenic_low_penetrance`
- `established_risk_allele`, `likely_risk_allele`, `uncertain_risk_allele`
- `affects`, `association`, `drug_response`, `confers_sensitivity`
- `protective`, `other`, `conflicting_interpretations_of_pathogenicity`, `not_provided`

### FAVOR Functional Annotations
| Node | Description |
|------|-------------|
| `annotation/info/FunctionalAnnotation/VarInfo` | Functional consequence |
| `annotation/info/FunctionalAnnotation/genecode_comprehensive_info` | Gene information |
| `annotation/info/FunctionalAnnotation/clnsig` | ClinVar clinical significance |
| `annotation/info/FunctionalAnnotation/clndn` | ClinVar disease name |
| `annotation/info/FunctionalAnnotation/bravo_af` | Bravo allele frequency |
| `annotation/info/FunctionalAnnotation/apc_protein_function_v3` | Protein function score |

## Usage

### DNAnexus Platform

#### Step 1: VCF Merge
```bash
dx run step1_vcf_merge \
  -ivcfs=sample1.vcf.gz \
  -ivcfs=sample2.vcf.gz \
  -o merged_output
```

#### Step 2: VCF to GDS Conversion
```bash
dx run vcf2gds \
  -ivcf_file=merged_chr1.vcf.gz \
  -igds_filename=merged_chr1.gds
```

#### Step 3: Functional Annotation
```bash
dx run favorannotator \
  -igds=merged_chr1.gds \
  -ichromosome=1 \
  -iuse_compression=TRUE
```

#### Step 4: Variant Prioritization
```bash
dx run step4_impact_prioritization \
  -igenelist=GeneList.txt \
  -igds_files=favor_merged_chr1.gds \
  -igds_files=favor_merged_chr2.gds
```

### Local Execution

For Step 4 local execution:
```bash
cd step4_impact_prioritization/resources/home/dnanexus/
Rscript IMPACT-prioritization.r --genelist GeneList.txt --outprefix anno_merged_
```

## Integration with IMPACT-VIS

The output `*_SNV_IMPACT.gds` files are designed for direct use with [IMPACT-VIS](https://boehlernick.github.io/IMPACT-VIS), an interactive R Shiny application for variant visualization.

### Preparing Data for IMPACT-VIS
1. Place output files in the IMPACT-VIS `app/data/` directory
2. Follow the naming convention: `{sample_id}_SNV_IMPACT.gds`
3. See the [IMPACT-VIS Data Preparation Guide](https://boehlernick.github.io/IMPACT-VIS/guides/02-data-preparation.html) for details

### IMPACT-VIS Capabilities
- Interactive visualization of prioritized variants
- Filtering by tier, ClinVar significance, allele frequency
- Gene-based and genomic region filtering
- Publication-ready plots with Plotly
- Persistent annotation states for sample curation

## Tiering System

Variants are assigned to tiers based on evidence strength:

| Tier | Criteria | Score Formula |
|------|----------|---------------|
| **Tier 1** | Pathogenic/Likely Pathogenic in ClinVar + gene match | 80 + 20 × globalScore |
| **Tier 2** | Frameshift, stopgain mutations | 60 + 40 × globalScore |
| **Tier 3** | Nonsynonymous, nonframeshift, stoploss | 20 + 80 × globalScore |
| **Tier 4** | APC protein function evidence | 100 × (0.5 × APC + 0.5 × globalScore) |

## Gene-Disease Association File

The `GeneList.txt` file should contain phenotype-specific gene associations from [Open Targets](https://www.opentargets.org/):

```
symbol	globalScore
GJB2	0.858985237
OTOF	0.850795742
MYO6	0.845566359
...
```

Generate this file by:
1. Querying Open Targets for your phenotype of interest
2. Exporting gene associations with global scores
3. Formatting as tab-separated with header row

## Refactor Planning and Contracts

The IMPACT-SNV pipeline is undergoing a phased refactor to replace the legacy FAVORannotator-based annotation step with a modern FAVOR-CLI-backed adapter, while preserving the existing prioritization logic and output contract.

For developers and contributors, comprehensive documentation is available:

- **[Pipeline Contract](docs/pipeline_contract.md)** - Current IMPACT-SNV workflow, inputs, outputs, and step-specific behavior
- **[FAVOR-CLI Refactor Plan](docs/favorcli_refactor_plan.md)** - Staged modernization strategy and architecture
- **[Annotation Compatibility Contract](docs/annotation_compatibility_contract.md)** - Required annotation fields and compatibility schema for Step 4
- **[Output GDS Contract](docs/output_gds_contract.md)** - Final output file structure and naming requirements
- **[Chromosome Handling Contract](docs/chromosome_handling_contract.md)** - Current (chr1-22) and target (chr1-22, X, Y) chromosome support
- **[FAVOR-CLI Schema Discovery Protocol](docs/favorcli_schema_discovery.md)** - Guidelines for empirical discovery of FAVOR-CLI output schema
- **[Documentation Index](docs/README.md)** - Overview of all refactor and contract documentation

These documents are designed to support code review, automated agents, and maintainability during the refactor. See [.github/copilot-instructions.md](.github/copilot-instructions.md) for additional GitHub Copilot coding agent instructions specific to this repository.

## References

### Tools and Databases
- [SeqArray](https://bioconductor.org/packages/SeqArray/) - Efficient storage of sequence data
- [FAVOR](https://favor.genohub.org/) - Functional Annotation of Variants Online Resource
- [Open Targets](https://www.opentargets.org/) - Gene-disease association platform
- [DNAnexus](https://www.dnanexus.com/) - Cloud-based genomics platform

### Related IMPACT Modules
- [IMPACT-VIS](https://github.com/boehlernick/IMPACT-VIS) - Interactive visualization
- [IMPACT-SV](https://github.com/boehlernick/IMPACT-SV) - Structural variant processing
- [IMPACT-CNV](https://github.com/boehlernick/IMPACT-CNV) - Copy number variant processing

## Citation

If you use this pipeline, please cite:

```bibtex
@article{impact-TBA,
  title={Integrated Mapping of Phenotype-Associated Candidate Targets for 
         interpretation and prioritization of genomic variants},
  authors={Boehler, N. and Cheng, H. Y. M.},
  journal={TBA},
  year={TBA}
}
```

### FAVOR Citation
> Zhou H., et al. (2023). FAVOR: functional annotation of variants online resource and annotator for variation across the human genome. *Nucleic Acids Research*, 51(D1), D1300-D1311. [DOI: 10.1093/nar/gkac966](https://doi.org/10.1093/nar/gkac966)

## License

This project is licensed under the terms specified in the [LICENSE](LICENSE) file.

## Support

- **Issues**: Report bugs or request features on the [GitHub Issues](https://github.com/boehlernick/IMPACT-SNV/issues) page
- **Documentation**: See individual step READMEs for detailed usage
- **IMPACT-VIS Help**: Visit the [IMPACT-VIS documentation](https://boehlernick.github.io/IMPACT-VIS)

