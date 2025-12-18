# Step 4: IMPACT Prioritization

This step performs variant prioritization on annotated GDS files using gene-disease association data and functional scores, producing output files compatible with [IMPACT-VIS](https://boehlernick.github.io/IMPACT-VIS) for interactive visualization.

## Features

- **Scoring**: Assigns pathogenicity scores (0-100) based on multiple evidence sources
- **Tiered Prioritization**: Classifies variants into Tier 1-4 based on evidence strength
- **ClinVar Integration**: Creates boolean flags for all ClinVar clinical significance categories
- **Per-Sample Output**: Generates individual `*_SNV_IMPACT.gds` files per sample
- **IMPACT-VIS Compatible**: Output format matches IMPACT-VIS input requirements

## Input Requirements

### Annotated GDS Files
- FAVOR-annotated GDS files from Step 3 (e.g., `merged_chr*.gds`)
- Must contain `annotation/info/FunctionalAnnotation/` nodes

### Gene-Disease Association File (`GeneList.txt`)
Tab-separated file with phenotype-specific gene associations:

```
symbol	globalScore
GJB2	0.858985237
OTOF	0.850795742
MYO6	0.845566359
```

**Columns:**
- `symbol`: Gene symbol (HGNC)
- `globalScore`: Open Targets association score (0-1)

Generate from [Open Targets](https://www.opentargets.org/) for your phenotype of interest. Export the data as .TSV files and save the first two columns as `GeneList.txt`.

## Output Format

### Output Files
- `{sample_id}_SNV_IMPACT.gds` - One file per sample

### GDS Annotations Added

#### IMPACT Scores
| Node | Type | Description |
|------|------|-------------|
| `annotation/info/impact_score` | numeric | Pathogenicity score (0-100) |
| `annotation/info/impact_score_calc` | character | Tier and calculation method |
| `annotation/info/tier` | integer | Priority tier (1-4) |

#### ClinVar Flags
Boolean indicators under `annotation/info/clnsig_flags/`:

| Flag | Description |
|------|-------------|
| `pathogenic` | ClinVar Pathogenic |
| `likely_pathogenic` | ClinVar Likely Pathogenic |
| `uncertain_significance` | ClinVar VUS |
| `likely_benign` | ClinVar Likely Benign |
| `benign` | ClinVar Benign |
| `pathogenic_low_penetrance` | Pathogenic with low penetrance |
| `likely_pathogenic_low_penetrance` | Likely Pathogenic with low penetrance |
| `established_risk_allele` | Established risk allele |
| `likely_risk_allele` | Likely risk allele |
| `uncertain_risk_allele` | Uncertain risk allele |
| `affects` | Affects phenotype |
| `association` | Disease association |
| `drug_response` | Drug response |
| `confers_sensitivity` | Confers sensitivity |
| `protective` | Protective variant |
| `conflicting_interpretations_of_pathogenicity` | Conflicting interpretations |
| `other` | Other/unrecognized category |
| `not_provided` | No ClinVar annotation |

## Tiering System

| Tier | Criteria | Score Formula |
|------|----------|---------------|
| **1** | ClinVar Pathogenic/Likely Pathogenic + gene match | 80 + 20 × globalScore |
| **2** | Frameshift insertion/deletion, stopgain | 60 + 40 × globalScore |
| **3** | Nonsynonymous SNV, nonframeshift indel, stoploss | 20 + 80 × globalScore |
| **4** | APC protein function score > 1 | 100 × (0.5 × APC + 0.5 × globalScore) |

## Scripts

| File | Description |
|------|-------------|
| `IMPACT-prioritization.r` | Main prioritization script |
| `seqarray_append.R` | Post-processing: adds tier, clnsig flags, renames score fields |

## Usage

### Local Execution

```bash
cd step4_impact_prioritization/resources/home/dnanexus/

# Run with default settings
Rscript IMPACT-prioritization.r

# Run with custom options
Rscript IMPACT-prioritization.r \
  --genelist GeneList.txt \
  --outprefix anno_merged_ \
  --prefix merged_chr \
  --pattern "merged_chr(.*)\\.gds"
```

**Options:**
| Option | Default | Description |
|--------|---------|-------------|
| `--genelist`, `-g` | `GeneList.txt` | Gene-disease association file |
| `--outprefix`, `-o` | `anno_merged_` | Output file prefix |
| `--prefix`, `-p` | `merged_chr` | Input GDS file prefix |
| `--pattern` | `merged_chr(.*)\\.gds` | Regex to extract chromosome |

### DNAnexus Platform

```bash
dx run step4_impact_prioritization \
  -igenelist=GeneList.txt \
  -igds_files=favor_merged_chr1.gds \
  -igds_files=favor_merged_chr2.gds \
  -igds_files=favor_merged_chr3.gds
```

**Applet Inputs:**
| Input | Type | Required | Description |
|-------|------|----------|-------------|
| `genelist` | file | Yes | Gene-disease association file (`.txt`) |
| `gds_files` | array:file | Yes | Annotated GDS files (`.gds`) |

**Applet Outputs:**
| Output | Type | Description |
|--------|------|-------------|
| `snv_impact_gds` | array:file | Per-sample `*_SNV_IMPACT.gds` files |
| `log_file` | file | Processing log |

## Integration with IMPACT-VIS

The output files are designed for direct use with IMPACT-VIS:

1. Copy `*_SNV_IMPACT.gds` files to `IMPACT-VIS/app/data/`
2. Launch IMPACT-VIS
3. Select samples from the dropdown

See the [IMPACT-VIS Data Preparation Guide](https://boehlernick.github.io/IMPACT-VIS/guides/02-data-preparation.html) for details.

## Dependencies

### R Packages (CRAN)
- `optparse`, `rlang`, `cli`, `dplyr`, `stringr`, `stringi`
- `parallel`, `readr`, `digest`, `tidyr`

### R Packages (Bioconductor)
- `SeqArray`, `SeqVarTools`, `gdsfmt`

## See Also

- [Main IMPACT-SNV README](../README.md) - Full pipeline overview
- [IMPACT-VIS Documentation](https://boehlernick.github.io/IMPACT-VIS) - Visualization platform
- [Open Targets](https://www.opentargets.org/) - Gene-disease associations
