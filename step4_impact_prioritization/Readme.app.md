# IMPACT-SNV Step 4: Variant Prioritization DNAnexus Applet

This applet runs the IMPACT-prioritization module to score and prioritize variants in annotated GDS files using gene-disease association data and functional scores.

## Inputs
- **Gene-disease association file** (`genelist`): Tab-delimited file (e.g., `GeneList.txt`)
- **Annotated GDS files** (`gds_files`): One or more annotated GDS files (e.g., `merged_chr*.gds`)

## Outputs
- All `<sample_id>_SNV_IMPACT.gds` files (one per sample)

## Usage Example
```
dx run step4_impact_prioritization -igenelist=GeneList.txt -igds_files=merged_chr1.gds, merged_chr2.gds,...
```

## Notes
- The applet will process all provided GDS files and output merged, prioritized GDS files for each sample.
- See the main README for more details on the scoring and tiering logic.
