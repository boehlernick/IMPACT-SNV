# IMPACT-SNV Release Readiness Check (2026-05-30)

## Scope

This check reflects the current release-candidate baseline, not a final `1.0.0` release.

Current package version under test: `1.0.0rc1`

Supported release path under this check:

- FAVOR annotation-only inputs
- `impact-snv extract-genotypes`
- `impact-snv build-gds`
- `impact-snv finalize-gds`
- `impact-snv validate-gds`

## Executed Validation Matrix

### 1. Focused CLI end-to-end regression

Command:

```bash
source .venv/bin/activate
pytest -q tests/test_release_cli_e2e.py
```

Result: **PASS** (`1 passed in 171.41s`)

Coverage highlights:
- release CLI path exercised through `build-gds`, `finalize-gds`, `validate-gds`, and `qc-build`
- fixture scope uses `tests/genebreaker_vcfs/merged_vcf_test`
- validates the package-owned path rather than the legacy DNAnexus wrapper flow

---

### 2. Packaging build and artifact content checks

Commands:

```bash
source .venv/bin/activate
python -m build
python -m zipfile -l dist/impact_snv-1.0.0rc1-py3-none-any.whl | grep 'impact_snv/resources/'
tar -tf dist/impact_snv-1.0.0rc1.tar.gz | grep 'impact_snv/resources/'
```

Result: **PASS**

Artifacts built:
- `dist/impact_snv-1.0.0rc1.tar.gz`
- `dist/impact_snv-1.0.0rc1-py3-none-any.whl`

Packaged resource checks: **PASS**
- `impact_snv/resources/favor_flat_to_seqarray_gds.R`
- `impact_snv/resources/impact_prioritize_gds.R`
- `impact_snv/resources/add_impact_vis_compat_nodes.R`

---

### 3. Fresh all-samples GeneBreaker build smoke

Evidence:

```bash
tests/genebreaker_vcfs/merged_vcf_test/gds_build_from_extracted_genotypes/build_manifest.json
```

Result: **PASS with expected warnings**

Observed outcomes:
- build manifest version: `1.0.0rc1`
- six samples built: `Case1_proband`, `Case1_mother`, `Case1_father`, `Case2_proband`, `Case2_mother`, `Case2_father`
- chromosome coverage recorded in manifest: `1-22`, `X`, `Y`
- six merged pre-prioritization GDS outputs were generated successfully
- `flat_vs_merged_delta = 0` for all samples
- QC warnings: `14`, all from expected zero-row autosome slices in the Case2 fixture context

---

### 4. Fresh all-samples finalize and final validation

Command:

```bash
source .venv/bin/activate
impact-snv finalize-gds --input-dir tests/genebreaker_vcfs/merged_vcf_test/gds_build_from_extracted_genotypes/gds_merged --gene-list tests/genebreaker_vcfs/merged_vcf_test/GeneList.txt --out-dir tests/genebreaker_vcfs/merged_vcf_test/final_from_extracted_genotypes --qc-mode warn --force
impact-snv validate-gds --input-dir tests/genebreaker_vcfs/merged_vcf_test/final_from_extracted_genotypes --qc-mode warn
```

Result: **PASS**

Key outcomes:
- finalize manifest written at `tests/genebreaker_vcfs/merged_vcf_test/final_from_extracted_genotypes/finalize_manifest.json`
- finalize manifest summary: `ok_count = 6`, `failed_count = 0`
- final validation returned `[OK]` for all six outputs
- finalized variant counts:
	- `Case1_father`: `503091`
	- `Case1_mother`: `507547`
	- `Case1_proband`: `502764`
	- `Case2_father`: `272459`
	- `Case2_mother`: `507534`
	- `Case2_proband`: `305305`

## Readiness Summary

Overall result: **PASS for release candidate `1.0.0rc1`**

### Green
- Focused package CLI end-to-end regression
- Packaging artifact generation and packaged R resource inclusion
- Fresh all-samples GeneBreaker build smoke from extracted genotypes
- Fresh all-samples `finalize-gds` run
- Fresh all-samples `validate-gds` run on finalized outputs

### Expected warnings / non-blocking notes
- The GeneBreaker build smoke still reports `14` warnings from expected zero-row autosome slices in the Case2 fixture set under `--qc-mode warn`
- Local R execution may emit `renv` out-of-sync notices and fall back to embedded required columns; these did not block GDS generation in this run

## Recommended Next Actions

1. Keep the repository version at `1.0.0rc1` until the actual final release decision; do not describe the repo as already released at `1.0.0`.
2. Preserve `build_manifest.json` and `finalize_manifest.json` from the fresh GeneBreaker run as merge and release evidence.

## Notes

This report supersedes the earlier representative-only finalize spot checks with a fresh all-samples GeneBreaker build, finalize, and validate pass using the extracted-genotypes release path.
