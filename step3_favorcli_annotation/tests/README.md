# Step 3 FAVOR-CLI Scaffold Tests

This directory is reserved for future fixtures and smoke tests for the experimental FAVOR-CLI annotation scaffold.

Current validation is intentionally lightweight and can be performed with:

- `bash step3_favorcli_annotation/src/run_favorcli.sh` in dry-run mode against a local chromosome-level GDS or VCF fixture.
- `python3 -m py_compile step1_vcf_merge/src/IMPACT_SNV_indel.py` for the already-updated chromosome handling path.
- `Rscript step3_favorcli_annotation/src/extract_variant_identity.R --input <gds> --input-type gds --output <tsv> --reference-genome-build GRCh38` when a fixture is available.

Future tests should cover:

- dry-run execution without FAVOR-CLI installed
- chromosome normalization for `1-22`, `X`, and `Y`
- canonical key generation
- validation failure when required configuration is missing