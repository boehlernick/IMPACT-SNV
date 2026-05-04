# Step 3: FAVOR-CLI Annotation Skeleton

This directory contains an experimental scaffolding module for a future FAVOR-CLI-based replacement for the legacy FAVORannotator Step 3.

This is not a drop-in replacement for [Step3_favorannotator-rap/](../Step3_favorannotator-rap/). It does not inject annotations into GDS, it does not implement production FAVOR-CLI schema mappings, and it does not claim scientific equivalence with the legacy annotation step.

The intended future flow is:

```text
GDS → variant identity → FAVOR-CLI → normalized adapter → GDS injection → Step 4
```

## Current Scope

The scaffold provides an execution boundary for a future adapter layer:

```text
GDS or VCF input
→ variant identity extraction
→ FAVOR-CLI-ready intermediate input
→ optional FAVOR-CLI invocation
→ raw FAVOR-CLI output location
```

Current behavior:

- Accepts a chromosome-level GDS file or a chromosome-level VCF/VCF.GZ file.
- Extracts a canonical variant identity table.
- Supports chromosome normalization for `1-22`, `X`, and `Y`.
- Validates configuration before any annotation step runs.
- Supports dry-run mode that does not require FAVOR-CLI to be installed.
- Does not mutate the input GDS file.
- Does not inject annotations into GDS.
- Does not implement production FAVOR-CLI schema mappings.

## Canonical Variant Identity

The intermediate table uses the following columns:

```text
variant.id
chromosome
position
ref
alt
canonical_variant_key
```

The canonical key format is:

```text
chromosome-position-ref-alt
```

Chromosomes are normalized to the contract used by [`docs/chromosome_handling_contract.md`](../docs/chromosome_handling_contract.md): `1-22`, `X`, and `Y` without silent conversion of numeric sex chromosome aliases.

## Configuration

The scaffold makes the following configuration explicit:

- `favor_cli_path`
- `favor_database_path` / `favor_database_file` on DNAnexus
- `favor_database_version`
- `reference_genome_build`
- `threads`
- `memory_budget_gb`
- `output_dir`
- `dry_run`

On DNAnexus, the applet accepts a chromosome-level GDS file or VCF file plus the configuration fields above.

## Dry Run

Dry-run mode is the recommended first validation path. It:

1. Validates configuration.
2. Extracts the variant identity table.
3. Writes intermediate files to the configured output directory.
4. Does not require FAVOR-CLI.
5. Does not run annotation.
6. Does not mutate the input GDS.

### Example

```bash
cd step3_favorcli_annotation/src

export FAVORCLI_GDS_FILE=/path/to/chr1.gds
export FAVORCLI_FAVOR_CLI_PATH=favor
export FAVORCLI_DATABASE_FILE=/path/to/favor/database
export FAVORCLI_DATABASE_VERSION=unknown
export FAVORCLI_REFERENCE_BUILD=GRCh38
export FAVORCLI_THREADS=4
export FAVORCLI_MEMORY_BUDGET_GB=8
export FAVORCLI_OUTPUT_DIR=/tmp/favorcli_skeleton
export FAVORCLI_DRY_RUN=true

bash run_favorcli.sh
```

For a local run, the wrapper writes:

- a variant identity table
- a validation report
- a run manifest describing the raw FAVOR-CLI output location reserved for future use

If you want to exercise the VCF boundary instead, set `FAVORCLI_VCF_FILE` to a local `.vcf` or `.vcf.gz` path and leave `FAVORCLI_GDS_FILE` empty.

## FAVOR-CLI Invocation

The current scaffold does not implement production FAVOR-CLI invocation. If `dry_run=false`, it captures `favor --version` when available and then fails clearly with a placeholder error rather than guessing the final command or schema mapping.

This is intentional. The exact FAVOR-CLI invocation and output schema should be discovered and documented before production integration.

## Validation

Lightweight validation is available through:

- `src/validate_favorcli_config.sh`
- `src/validate_favorcli_output.R`
- `Rscript src/extract_variant_identity.R --help` style argument parsing

Where real GDS fixtures are unavailable, manual dry-run validation should confirm that the module produces the intermediate table without mutating the source GDS.

## Not Implemented Yet

This scaffold intentionally does not include:

- FAVOR-CLI schema mapping
- GDS annotation injection
- Step 4 compatibility normalization
- legacy FAVORannotator replacement
- any change to Step 4 scoring or final output naming

## Files

- `dxapp.json` - Experimental DNAnexus applet definition
- `src/code.sh` - DNAnexus entrypoint
- `src/run_favorcli.sh` - Run orchestration and dry-run handling
- `src/extract_variant_identity.R` - Read-only variant identity extraction
- `src/validate_favorcli_config.sh` - Configuration validation helpers
- `src/validate_favorcli_output.R` - Variant identity validation helpers
- `tests/` - Placeholder for future fixtures and smoke tests