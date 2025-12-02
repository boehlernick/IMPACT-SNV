
# Smoke Test for IMPACT-SNV Step 1

This smoke test validates the containerized Step 1 (VCF merge + normalization + per-contig split) on a tiny synthetic dataset.

## Purpose
The smoke test ensures that the container successfully:
1. Discovers and prepares input VCFs (bgzipping and indexing as needed)
2. Merges multi-sample VCFs
3. Normalizes multiallelic variants
4. Splits output into exactly 25 canonical contigs (Chr1–Chr22, ChrX, ChrY, ChrMT)
5. Creates empty header-only VCFs for missing contigs
6. Produces valid bgzipped outputs with proper indexes

## Prerequisites
- Docker installed and available in your PATH.
- The Step 1 Dockerfile and entrypoint.sh present in the project root.

## Files
- `smoke/data/`: Contains three small subset VCFs:
  - `NA12892.subset.sorted.vcf.gz`
  - `NA18591.subset.sorted.vcf.gz`
  - `NA18596.subset.sorted.vcf.gz`
- `smoke/out/`: output directory (created by the test).
- `assert.sh`: post-run assertions.
- `Makefile`: build and run targets.

## Usage

```bash
make test-step1
```

This will:
1. Build the Docker image (`impact-step1:latest`)
2. Run the container with deterministic inputs (`--input-vcfs`)
3. Validate outputs with `smoke/assert.sh`

Expected outputs under `smoke/out/`:
- `smoke.vcf.gz` + `smoke.vcf.gz.tbi` (merged)
- `smoke.normalized.vcf.gz` + `smoke.normalized.vcf.gz.tbi` (normalized)
- `smoke.normalized.Chr1.vcf.gz` through `smoke.normalized.ChrMT.vcf.gz` (25 contigs)
- Index files for each contig (`.tbi` or `.csi`)
- Total: 50+ files (27 VCFs + 25+ indexes)

## Cleaning
```bash
make clean
```
Removes `smoke/out/`.

## Smoke Test Assertions
The `assert.sh` script verifies:
- Merged VCF exists and is readable
- Normalized VCF exists and is readable
- At least Chr1 and ChrX contig splits exist
- Exactly 25 contig files are present
- All outputs are valid bgzipped VCFs (can be indexed/read by bcftools)
- Index files exist for all outputs

## Troubleshooting

**Test fails with "Missing merged VCF":**
- Check Docker image built successfully: `docker images | grep impact-step1`
- Review container logs: `docker run ... 2>&1 | tail -20`

**Test fails with "Missing per-contig splits":**
- Verify bcftools is installed in the container (check Dockerfile)
- Check that input VCFs contain the expected contigs

**Empty contig VCFs created:**
- Normal if the input data doesn't contain all 25 contigs
- The script creates header-only files for missing contigs to ensure downstream tools accept the output
