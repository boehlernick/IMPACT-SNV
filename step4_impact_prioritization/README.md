# Step 4: IMPACT Prioritization

This step performs variant prioritization on annotated GDS files using gene-disease association data and functional scores. The module has been substantially updated to:
- Support robust scoring and tiering of variants using gene-disease association and functional annotation.
- Output per-sample, per-chromosome GDS files and automatically merge them into final *SNV_IMPACT.gds files.
- Provide detailed logging and error handling.
- Be fully compatible as a DNAnexus applet for cloud-based execution.

## Input
- Annotated GDS files (e.g., `merged_chr*.gds`)
- Gene-disease association file: `GeneList.txt` (default, can be changed with --genelist)

## Main Script
- `IMPACT-prioritization.r`: R script to score and prioritize variants by pathogenicity and gene relevance.

## Usage (Local)
Run in a directory containing the required GDS files and `GeneList.txt`:

```sh
Rscript IMPACT-prioritization.r --genelist GeneList.txt --outprefix anno_merged_
```

- `--genelist` (optional): Path to gene-disease association file (default: `GeneList.txt`)
- `--outprefix` (optional): Prefix for output GDS files (default: `anno_merged_`)

The script will output new GDS files with prioritization scores and filtered variants for downstream analysis. Final merged files will be named `<sample_id>_SNV_IMPACT.gds`.

## DNAnexus Applet
This module is available as a DNAnexus applet (`step4_IMPACT_prioritization`).

### Applet Inputs
- **GeneList.txt**: Gene-disease association file (tab-delimited, required)
- **GDS files**: Annotated GDS files (e.g., `merged_chr*.gds`, required)

### Applet Outputs
- All `<sample_id>_SNV_IMPACT.gds` files (one per sample)

### Example DNAnexus Command
```
dx run step4_IMPACT_prioritization -igenelist=GeneList.txt -igds_files=merged_chr1.gds, merged_chr2.gds,...
```

See the main repository README for the full pipeline and context.
