# IMPACT-SNV

**Integrated Mapping of Phenotype-Associated Candidate Targets for SNV/Indel Analysis**

This repository contains the IMPACT-SNV pipeline, which processes and prioritizes single nucleotide variants (SNVs) and indels using FAVOR-based annotation and phenotype-specific gene-disease associations.

## Overview

IMPACT-SNV is part of the broader IMPACT framework for phenotype-configurable interpretation of genomic variants. This module specifically handles SNV/Indel processing and produces output files compatible with [IMPACT-VIS](https://boehlernick.github.io/IMPACT-VIS) for interactive visualization and analysis.

## Current Repo State

This repository currently contains both:
- the maintained `impact_snv/` Python package and `impact-snv` CLI for local and refactor-era execution
- the original DNAnexus step directories retained for backward compatibility and platform packaging

Implemented package CLI commands are:
- `sanitize-vcfs`
- `merge`
- `favor-ingest`
- `favor-annotate`
- `extract-genotypes`
- `build-gene-lists`
- `build-gds`
- `finalize-gds`
- `validate-gds`
- `qc-build`

The umbrella `impact-snv run` command is still planned and is not implemented yet. Current focused automated coverage includes contract tests plus a real CLI integration test for `build-gds`, `finalize-gds`, `validate-gds`, and `qc-build`.


### Why FAVOR?

IMPACT-SNV currently supports two Step 3 annotation backends during the refactor window:
- legacy FAVORannotator compatibility (`legacy-favorannotator` backend)
- FAVOR CLI (`favor-cli` backend)

FAVOR was selected because it provides a large precomputed and harmonized annotation resource suitable for WGS-scale variant prioritization. The annotation set includes population allele frequencies, transcript consequences, conservation metrics, ClinVar assertions, regulatory annotations, and aggregate protein impact metrics.

In IMPACT-SNV, FAVOR annotations are appended to GDS/aGDS files generated from VCF inputs, allowing large variant datasets to be stored, queried, and prioritized efficiently. This design fits the workflow’s goal of scalable phenotype-driven prioritization across WGS-derived variant calls.

A FAVOR-CLI-compatible path now exists in the repository through the backend-based `impact-snv favor-annotate` command and the RAP-oriented `step3_favorcli_annotation/` module. Legacy FAVORannotator support remains in place for backward compatibility while the FAVOR-CLI path is hardened against the existing Step 4 and output-GDS contracts.

For the supported v1.0.0 release path, FAVOR is currently treated as annotation-only. Native FAVOR `.cohort` or FAVOR-generated genotype outputs remain experimental/discovery-only and are not part of the release gate; the supported path is `favor-annotate` plus `impact-snv extract-genotypes`.

## Pipeline Architecture

The repository still follows the historical four-step IMPACT-SNV model, and the legacy step folders remain available as DNAnexus applets.

```
VCF Files → [Step 1: Merge] → [Step 2: VCF2GDS] → [Step 3: FAVOR Annotate] → [Step 4: Prioritize] → *_SNV_IMPACT.gds
```

The maintained package CLI exposes the same workflow through more granular local commands:

```
VCF Files → impact-snv sanitize-vcfs (when sample IDs collide) → impact-snv merge → impact-snv favor-ingest / favor-annotate + impact-snv extract-genotypes → impact-snv build-gene-lists → impact-snv build-gds → impact-snv finalize-gds → impact-snv validate-gds / qc-build
```

| Step | Folder | Description | Input | Output |
|------|--------|-------------|-------|--------|
| 1 | `step1_vcf_merge/` | Merges multiple VCF files into chromosome-separated files | `.vcf`, `.vcf.gz` | `merged_chr*.vcf.gz` |
| 2 | `step2_vcf2gds/` | Converts VCF to GDS format for efficient processing | `.vcf.gz` | `merged_chr*.gds` |
| 3 | `Step3_favorannotator-rap/`, `step3_favorcli_annotation/`, `impact_snv/favor/` | Annotates variants using the selected Step 3 backend | legacy `.gds`, FAVOR-ingested directories, or backend dry-run inputs | canonical `<out_prefix>.annotated` outputs plus compatibility-staged inputs |
| 4 | `step4_impact_prioritization/`, `impact_snv/gene_lists.py`, `impact_snv/gds/` | Builds phenotype-driven per-sample GeneLists, scores, finalizes, validates, and QC-checks per-sample outputs | annotated legacy GDS or pre-prioritization per-sample GDS plus a single `GeneList.txt` or a per-sample gene-list manifest | `*_SNV_IMPACT.gds`, manifests, and QC summaries |

## Requirements

### Python Package And Local CLI
- **Python** >= 3.10
- Python package dependencies declared in `pyproject.toml`: `pandas`, `pyarrow`
- **Rscript** with the packaged GDS path dependencies:
  - `optparse`, `arrow`, `SeqArray`, `gdsfmt`, `stringi`
- **BCFtools**, **bgzip**, and **tabix** when using `impact-snv merge` or `impact-snv extract-genotypes`
- A **FAVOR CLI** binary named `favor` when using the `favor-cli` backend

### Legacy Applet Extras
- `SeqVarTools`, `dplyr`, and `tidyr` remain required by legacy R/DNAnexus components
- **Rust** and **xsv** are only required for the legacy `Step3_favorannotator-rap/` path

### Platform
The repository still contains **DNAnexus Research Analysis Platform** applets under the step folders, each with its own `dxapp.json`. The packaged `impact-snv` CLI can also be run locally and is the maintained surface for the refactor-era build, finalize, validate, and QC workflow.

## Input Requirements

### Step 1: VCF Merge
- **Input**: One or more VCF or VCF.GZ files
- **Output**: Chromosome-separated VCF.GZ files (`merged_chr1.vcf.gz`, `merged_chr2.vcf.gz`, etc.)

### Step 2: VCF to GDS
- **Input**: VCF.GZ file from Step 1
- **Output**: GDS file (`merged_chr*.gds`)

### Step 3: FAVOR Annotation
- **`favor-cli` backend input**: FAVOR-ingested directory (`<out_prefix>.ingested`) produced by `impact-snv favor-ingest`. The input VCF must have a `.csi` or `.tbi` index alongside it.
- **`legacy-favorannotator` backend input**: existing legacy annotated/genotype directories for compatibility staging
- **`favor-cli-skeleton` backend input**: GDS or VCF for dry-run planning
- **Output**: Canonical `<out_prefix>.annotated` directory written to `<out_dir>`, consumed by `impact-snv build-gds`; genotype extraction is handled separately by `impact-snv extract-genotypes`

The package CLI `legacy-favorannotator` backend stages existing legacy outputs for downstream compatibility. It does not locally invoke the old DNAnexus FAVORannotator applet.

### Step 4: IMPACT Prioritization And Finalization
- **Input**:
  - Legacy annotated chromosome GDS files from Step 3, or pre-prioritization per-sample GDS files from `impact-snv build-gds`
  - Gene-disease association input, either:
    - one shared `GeneList.txt`, or
    - a per-sample gene-list manifest produced by `impact-snv build-gene-lists`
  - Each `GeneList.txt`-style file remains tab-separated with columns:
    - `symbol`: Gene symbol (e.g., `GJB2`, `OTOF`)
    - `globalScore`: Open Targets association score (0-1)
- **Output**: Per-sample final GDS files (`{sample_id}_SNV_IMPACT.gds`)

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
| `annotation/info/impact_score` | numeric | Prioritization score (0-100) |
| `annotation/info/impact_score_calc` | character | Tier and calculation formula |
| `annotation/info/tier` | integer | Priority tier (1-4) |
| `annotation/info/scoring_gene` | character | Gene selected for IMPACT scoring |
| `annotation/info/scoring_gene_score` | numeric | Gene-level score used in the tier formula |

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
| `annotation/info/FunctionalAnnotation/VarInfo` | IMPACT-VIS variant key in `chr-pos-ref-alt` format |
| `annotation/info/FunctionalAnnotation/genecode_comprehensive_info` | Gene information |
| `annotation/info/FunctionalAnnotation/clnsig` | ClinVar clinical significance |
| `annotation/info/FunctionalAnnotation/clndn` | ClinVar disease name |
| `annotation/info/FunctionalAnnotation/bravo_af` | Bravo allele frequency |
| `annotation/info/FunctionalAnnotation/aloft_prediction` | ALOFT compatibility field expected by downstream readers |
| `annotation/info/FunctionalAnnotation/apc_protein_function_v3` | Protein function score |

### Annotation Compatibility And Provenance Nodes
When building GDS from flat parquet, IMPACT-SNV now records explicit fallback/provenance metadata under:

- `annotation/info/IMPACT_AnnotationCompatibility/fallback_flags`
- `annotation/info/IMPACT_AnnotationCompatibility/fallback_count`
- `annotation/info/IMPACT_AnnotationProvenance/adapter_version`
- `annotation/info/IMPACT_AnnotationProvenance/compatibility_contract_version`

These nodes make defaulted/backfilled values auditable rather than implicit.

## Usage

### DNAnexus Applets

The original DNAnexus step implementations remain in:
- `step1_vcf_merge/`
- `step2_vcf2gds/`
- `Step3_favorannotator-rap/`
- `step3_favorcli_annotation/`
- `step4_impact_prioritization/`

Use each step directory's README and `dxapp.json` for platform-specific invocation details. The root README focuses on the current repository shape and the maintained local CLI surface.

### Local Execution

Install the package in your environment and inspect the available commands:

```bash
python -m pip install -e .
impact-snv --help
```

When input VCFs come from repeated trio-style fixtures with sample IDs like `proband`, `mother`, and `father`, sanitize them before merge so downstream GDS and IMPACT-VIS sample IDs stay unique and stable:

```bash
impact-snv sanitize-vcfs \
  --vcfs tests/production_test/*.vcf.gz \
  --out-dir out/sanitized_inputs

impact-snv merge \
  --vcf-manifest out/sanitized_inputs/sanitized_vcfs.tsv \
  --out-vcf out/merged.vcf.gz
```

By default the sanitizer derives a case prefix from the filename, so trio sample IDs like `proband`, `mother`, and `father` become names such as `Case1_proband`, `Case1_mother`, and `Case1_father`.

To build phenotype-specific per-sample GeneLists from Open Targets, provide a sample manifest and a phenotype manifest. The resulting `sample_gene_lists.tsv` can be passed directly into `build-gds` and `finalize-gds`:

```bash
impact-snv build-gene-lists \
  --samples-manifest tests/production_test/samples.csv \
  --phenotypes-manifest tests/production_test/phenotypes.csv \
  --out-dir out/gene_lists
```

Step 3 runs in two stages. First, ingest the merged VCF into FAVOR's variant-set format (`<out_prefix>.ingested`). The input VCF must have a `.csi` or `.tbi` index alongside it:

```bash
impact-snv favor-ingest \
  --input-vcf out/merged.vcf.gz \
  --out-dir out/ \
  --out-prefix case1
```

Then annotate the ingested variant set against the FAVOR database. This is the long-running step — expect several hours for WGS-scale cohorts. Future releases will seek to parallelize this operation to increase efficiency. The default backend is `favor-cli`; `--backend` may be omitted:

```bash
impact-snv favor-annotate \
  --ingested-dir out/case1.ingested \
  --out-dir out/ \
  --out-prefix case1
```

Legacy compatibility staging (for pre-existing FAVORannotator outputs):

```bash
impact-snv favor-annotate \
  --backend legacy-favorannotator \
  --legacy-annotated-dir path/to/legacy.annotated \
  --legacy-genotypes-dir path/to/legacy.genotypes \
  --legacy-stage-mode symlink \
  --out-dir out/ \
  --out-prefix case1
```

This backend stages pre-existing legacy outputs only. It does not execute the DNAnexus `Step3_favorannotator-rap/` applet locally.

`impact-snv favor-annotate` does not replace genotype extraction. For the maintained local release path, extract genotypes separately and then build, finalize, validate, and QC the per-sample GDS outputs:

```bash
impact-snv extract-genotypes \
  --input-vcf out/merged.vcf.gz \
  --out-dir out/genotypes

impact-snv build-gds \
  --annotated-dir out/merged.annotated \
  --genotypes-dir out/genotypes \
  --gene-list-manifest out/gene_lists/sample_gene_lists.tsv \
  --out-dir out/build \
  --all-samples

impact-snv finalize-gds \
  --input-dir out/build/gds_merged \
  --gene-list-manifest out/gene_lists/sample_gene_lists.tsv \
  --out-dir out/final

impact-snv validate-gds \
  --input-dir out/final \
  --qc-mode strict

impact-snv qc-build \
  --manifest out/build/build_manifest.json \
  --out-dir out/qc \
  --qc-mode strict
```

For a concrete tested example of the build/finalize/validate/QC path, see `tests/test_release_cli_e2e.py`.

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

## Gene-Disease Association Files

Each `GeneList.txt`-style file should contain phenotype-specific gene associations from [Open Targets](https://www.opentargets.org/):

```
symbol	globalScore
GJB2	0.858985237
OTOF	0.850795742
MYO6	0.845566359
...
```

The maintained CLI can now generate these files directly from per-sample phenotype manifests:
1. `impact-snv build-gene-lists` reads a samples manifest and a phenotype manifest.
2. For each sample, it queries Open Targets for every phenotype or HPO term assigned to that sample.
3. It retains the union of all returned gene symbols and keeps the highest `globalScore` observed for each symbol across that sample's phenotype queries.
4. It writes one `GeneList.txt`-compatible file per sample plus a `sample_gene_lists.tsv` manifest that `build-gds` and `finalize-gds` can consume.

When phenotype text is provided without an ontology ID, IMPACT-SNV resolves it through the Open Targets search endpoint, prefers an exact normalized name match when available, and otherwise falls back to the top disease/phenotype search hit. The selected disease IDs and match strategies are recorded in the JSON manifest for auditability.

## Refactor Planning and Contracts

The repository is in a phased, backwards-compatible refactor. Legacy FAVORannotator-era applets are still present, while the packaged `impact-snv` CLI now owns the maintained local interfaces for Step 3 backend selection, GDS build/finalize/validate, and build QC. The goal remains to preserve the existing prioritization logic and output contract while expanding support for FAVOR-CLI-backed annotation.

For developers and contributors, comprehensive documentation is available:

### Development Guidelines
- **[Copilot Instructions](.github/copilot-instructions.md)** - Development guardrails, FAVOR-CLI principles, and critical guardrails for coding agents

### Architectural References
- **[FAVOR-CLI Refactor Plan](docs/favorcli_refactor_plan.md)** - Staged modernization strategy (10-milestone roadmap) and overall architecture
- **[Pipeline Contract](docs/pipeline_contract.md)** - Current IMPACT-SNV workflow, inputs, outputs, and step-specific behavior

### Contract Specifications
- **[Annotation Compatibility Contract](docs/annotation_compatibility_contract.md)** - Required annotation fields and compatibility schema for Step 4
- **[Output GDS Contract](docs/output_gds_contract.md)** - Final output file structure and naming requirements
- **[Chromosome Handling Contract](docs/chromosome_handling_contract.md)** - Current chromosome normalization and supported `1-22`, `X`, `Y` handling

### Discovery and Validation
- **[FAVOR-CLI Schema Discovery Protocol](docs/favorcli_schema_discovery.md)** - Guidelines for empirical discovery of FAVOR-CLI output schema
- **[Release Readiness Check](docs/release_readiness_2026-05-30.md)** - Current regression, packaging, finalize, and QC validation summary

### Documentation Index
- **[Complete Documentation Index](docs/README.md)** - Overview of all refactor and contract documentation

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

