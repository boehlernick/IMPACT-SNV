# Step 4: IMPACT Prioritization

This step performs variant prioritization on annotated GDS files using gene-disease association data and functional scores.

## Input
- Annotated GDS files (e.g., `favor_merged_chr*.gds`)
- Gene-disease association file: `GeneList.txt` (default, can be changed with --gda)

## Main Script
- `IMPACT-prioritization.r`: R script to score and prioritize variants by pathogenicity and gene relevance.

## Usage
Run in a directory containing the required GDS files and `GeneList.txt`:

```sh
Rscript IMPACT-prioritization.r --gda GeneList.txt --outprefix anno_merged_
```

- `--gda` (optional): Path to gene-disease association file (default: `GeneList.txt`)
- `--outprefix` (optional): Prefix for output GDS files (default: `anno_merged_`)

The script will output new GDS files with prioritization scores and filtered variants for downstream analysis.

See the main repository README for the full pipeline and context.
