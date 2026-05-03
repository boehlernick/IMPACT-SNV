#!/bin/bash
set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

main() {
  local working_dir="/home/dnanexus"
  local output_dir="${output_dir:-out/favorcli_skeleton}"
  local gds_local=""
  local vcf_local=""
  local favor_database_local=""

  mkdir -p "$working_dir" "$output_dir"

  if [[ -n "${gds_file:-}" ]]; then
    dx download "$gds_file" -o "$working_dir/input.gds"
    gds_local="$working_dir/input.gds"
  fi

  if [[ -n "${vcf_file:-}" ]]; then
    dx download "$vcf_file" -o "$working_dir/input.vcf.gz"
    vcf_local="$working_dir/input.vcf.gz"
  fi

  if [[ -n "${favor_database_file:-}" ]]; then
    dx download "$favor_database_file" -o "$working_dir/favor_database_file"
    favor_database_local="$working_dir/favor_database_file"
  fi

  export FAVORCLI_GDS_FILE="$gds_local"
  export FAVORCLI_VCF_FILE="$vcf_local"
  export FAVORCLI_DATABASE_FILE="$favor_database_local"
  export FAVORCLI_OUTPUT_DIR="$output_dir"
  export FAVORCLI_FAVOR_CLI_PATH="${favor_cli_path:-favor}"
  export FAVORCLI_DATABASE_VERSION="${favor_database_version:-}"
  export FAVORCLI_REFERENCE_BUILD="${reference_genome_build:-GRCh38}"
  export FAVORCLI_THREADS="${threads:-4}"
  export FAVORCLI_MEMORY_BUDGET_GB="${memory_budget_gb:-8}"
  export FAVORCLI_DRY_RUN="${dry_run:-true}"

  source "$SCRIPT_DIR/validate_favorcli_config.sh"
  validate_favorcli_config

  bash "$SCRIPT_DIR/run_favorcli.sh"

  mkdir -p out/variant_identity_table out/validation_report out/run_manifest
  cp "$output_dir/variant_identity.tsv" out/variant_identity_table/
  cp "$output_dir/favorcli_validation_report.txt" out/validation_report/
  cp "$output_dir/favorcli_run_manifest.json" out/run_manifest/

  dx-upload-all-outputs
}

main "$@"