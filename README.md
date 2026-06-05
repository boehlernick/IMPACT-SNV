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

FAVOR is used as an annotation backend for WGS-scale prioritization. In the maintained local path, annotations are produced first and genotypes are extracted separately.

The repository includes a hardened FAVOR-CLI path (`impact-snv favor-annotate`) while keeping legacy FAVORannotator compatibility for backward support. For the current release path, FAVOR `.cohort` or FAVOR-generated genotype outputs are not part of the supported release gate.

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

## Installation

### Prerequisites

Minimum required tools for the maintained local CLI path:

- **Python** >= 3.10
- **R** with `Rscript`
- **BCFtools**, **bgzip**, and **tabix** (required for `impact-snv merge` and `impact-snv extract-genotypes`)
- **FAVOR CLI** executable named `favor` on `PATH` when using backend `favor-cli`

Python package dependencies are declared in `pyproject.toml` (`pandas`, `pyarrow`).

R-side GDS scripts require packages including `optparse`, `arrow`, `SeqArray`, `gdsfmt`, and `stringi`.

### Install The Python CLI

From repository root:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -e .
impact-snv --help
```

### Install / Repair R Dependencies (renv)

This repository includes `renv.lock` and activates `renv` in packaged R scripts.

Recommended setup:

```bash
R -q -e "renv::restore(prompt = FALSE)"
```

If your cache contains broken symlinks or missing packages, repair then restore:

```bash
R -q -e "Sys.setenv(RENV_CONFIG_CACHE_SYMLINKS='FALSE'); renv::repair(); renv::restore(prompt = FALSE)"
```

Optional dependency sanity check:

```bash
R -q -e "library(optparse); library(arrow); library(SeqArray); library(gdsfmt)"
```

### Platform Notes

The repository still contains **DNAnexus Research Analysis Platform** applets under the step folders (each with its own `dxapp.json`).
The packaged `impact-snv` CLI is the maintained local interface for build/finalize/validate/QC workflows.

## Inputs And Outputs (Concise)

Maintained local workflow inputs:

- Variant input: one or more VCF/VCF.GZ files
- Step 3 annotation input: FAVOR-ingested directory (`<out_prefix>.ingested`) for backend `favor-cli`
- Genotype input for build stage: output directory from `impact-snv extract-genotypes`
- Gene scoring input: one shared `GeneList.txt` or a per-sample manifest from `impact-snv build-gene-lists`

Maintained local workflow outputs:

- Step 3: canonical `<out_prefix>.annotated` directory
- Build stage: per-sample `*_SNV_IMPACT.preprioritization.gds`
- Final stage: per-sample `{sample_id}_SNV_IMPACT.gds`

For full schema and contract details, see:

- `docs/pipeline_contract.md`
- `docs/annotation_compatibility_contract.md`
- `docs/output_gds_contract.md`
- `docs/chromosome_handling_contract.md`

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

#### Quickstart (Maintained Local Path)

This is the current supported local sequence:

```bash
impact-snv sanitize-vcfs \
  --vcfs tests/production_test/*.vcf.gz \
  --out-dir out/sanitized_inputs

impact-snv merge \
  --vcf-manifest out/sanitized_inputs/sanitized_vcfs.tsv \
  --out-vcf out/merged.vcf.gz

impact-snv favor-ingest \
  --input-vcf out/merged.vcf.gz \
  --out-dir out/ \
  --out-prefix case1

impact-snv favor-annotate \
  --ingested-dir out/case1.ingested \
  --out-dir out/ \
  --out-prefix case1

impact-snv extract-genotypes \
  --input-vcf out/merged.vcf.gz \
  --out-dir out/genotypes

impact-snv build-gene-lists \
  --samples-manifest tests/production_test/samples.csv \
  --phenotypes-manifest tests/production_test/phenotypes.csv \
  --out-dir out/gene_lists

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

#### Notes For Resuming Interrupted `build-gds` Runs

- `build-gds` skips existing per-chromosome flat parquet and GDS outputs unless `--force` is provided.
- Re-running the same command is safe for resume behavior.
- Use `--force` only when you want to overwrite generated outputs.

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

