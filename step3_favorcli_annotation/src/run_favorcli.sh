#!/bin/bash
set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

normalize_bool() {
  local value="${1:-false}"
  case "${value,,}" in
    1|true|yes|y)
      printf 'true'
      ;;
    *)
      printf 'false'
      ;;
  esac
}

write_manifest() {
  local manifest_path="$1"
  cat > "$manifest_path" <<EOF
{
  "dry_run": "$FAVORCLI_DRY_RUN",
  "input_gds": "$FAVORCLI_GDS_FILE",
  "input_vcf": "$FAVORCLI_VCF_FILE",
  "favor_cli_path": "$FAVORCLI_FAVOR_CLI_PATH",
  "favor_database_file": "$FAVORCLI_DATABASE_FILE",
  "favor_database_version": "$FAVORCLI_DATABASE_VERSION",
  "reference_genome_build": "$FAVORCLI_REFERENCE_BUILD",
  "threads": "$FAVORCLI_THREADS",
  "memory_budget_gb": "$FAVORCLI_MEMORY_BUDGET_GB",
  "raw_output_dir": "$FAVORCLI_OUTPUT_DIR/raw_favorcli_output",
  "status": "${2:-dry_run_completed}"
}
EOF
}

main() {
  local output_dir="$FAVORCLI_OUTPUT_DIR"
  local variant_identity_tsv="$output_dir/variant_identity.tsv"
  local validation_report="$output_dir/favorcli_validation_report.txt"
  local manifest_path="$output_dir/favorcli_run_manifest.json"
  local raw_output_dir="$output_dir/raw_favorcli_output"
  local dry_run

  dry_run=$(normalize_bool "$FAVORCLI_DRY_RUN")
  mkdir -p "$output_dir" "$raw_output_dir"

  if [[ -n "$FAVORCLI_GDS_FILE" ]]; then
    Rscript "$SCRIPT_DIR/extract_variant_identity.R" \
      --input "$FAVORCLI_GDS_FILE" \
      --input-type gds \
      --output "$variant_identity_tsv" \
      --reference-genome-build "$FAVORCLI_REFERENCE_BUILD"
  elif [[ -n "$FAVORCLI_VCF_FILE" ]]; then
    Rscript "$SCRIPT_DIR/extract_variant_identity.R" \
      --input "$FAVORCLI_VCF_FILE" \
      --input-type vcf \
      --output "$variant_identity_tsv" \
      --reference-genome-build "$FAVORCLI_REFERENCE_BUILD"
  else
    echo "Error: expected either a GDS file or a VCF file." >&2
    exit 1
  fi

  Rscript "$SCRIPT_DIR/validate_favorcli_output.R" \
    --input "$variant_identity_tsv" \
    --output "$validation_report" \
    --reference-genome-build "$FAVORCLI_REFERENCE_BUILD" \
    --dry-run "$dry_run"

  if [[ "$dry_run" == "true" ]]; then
    write_manifest "$manifest_path" "dry_run_completed"
    printf 'Dry-run complete. Variant identity table written to %s\n' "$variant_identity_tsv"
    printf 'Raw FAVOR-CLI output location reserved at %s\n' "$raw_output_dir"
    return 0
  fi

  if command -v "$FAVORCLI_FAVOR_CLI_PATH" >/dev/null 2>&1; then
    "$FAVORCLI_FAVOR_CLI_PATH" --version | tee "$output_dir/favorcli_version.txt"
  else
    echo "Error: FAVOR-CLI executable not found at '$FAVORCLI_FAVOR_CLI_PATH'." >&2
    write_manifest "$manifest_path" "missing_favor_cli"
    exit 1
  fi

  write_manifest "$manifest_path" "placeholder_not_implemented"
  echo "Error: FAVOR-CLI invocation is intentionally not implemented in this scaffold." >&2
  echo "Raw output would be written to: $raw_output_dir" >&2
  exit 1
}

main "$@"