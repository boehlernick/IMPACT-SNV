#!/bin/bash
set -euo pipefail

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

validate_favorcli_config() {
  local dry_run
  dry_run=$(normalize_bool "${FAVORCLI_DRY_RUN:-false}")

  if [[ -z "${FAVORCLI_GDS_FILE:-}" && -z "${FAVORCLI_VCF_FILE:-}" ]]; then
    echo "Error: provide either a GDS input or a VCF input." >&2
    exit 1
  fi

  if [[ -n "${FAVORCLI_GDS_FILE:-}" && ! -f "${FAVORCLI_GDS_FILE}" ]]; then
    echo "Error: GDS input not found: ${FAVORCLI_GDS_FILE}" >&2
    exit 1
  fi

  if [[ -n "${FAVORCLI_VCF_FILE:-}" && ! -f "${FAVORCLI_VCF_FILE}" ]]; then
    echo "Error: VCF input not found: ${FAVORCLI_VCF_FILE}" >&2
    exit 1
  fi

  if [[ -n "${FAVORCLI_DATABASE_FILE:-}" && ! -f "${FAVORCLI_DATABASE_FILE}" ]]; then
    echo "Error: FAVOR database file not found: ${FAVORCLI_DATABASE_FILE}" >&2
    exit 1
  fi

  if [[ -z "${FAVORCLI_OUTPUT_DIR:-}" ]]; then
    echo "Error: output directory is required." >&2
    exit 1
  fi

  if [[ -z "${FAVORCLI_REFERENCE_BUILD:-}" ]]; then
    echo "Error: reference genome build is required." >&2
    exit 1
  fi

  if ! [[ "${FAVORCLI_THREADS:-4}" =~ ^[0-9]+$ ]] || [[ "${FAVORCLI_THREADS:-4}" -lt 1 ]]; then
    echo "Error: threads must be a positive integer." >&2
    exit 1
  fi

  if ! [[ "${FAVORCLI_MEMORY_BUDGET_GB:-8}" =~ ^[0-9]+$ ]] || [[ "${FAVORCLI_MEMORY_BUDGET_GB:-8}" -lt 1 ]]; then
    echo "Error: memory budget must be a positive integer." >&2
    exit 1
  fi

  if [[ "$dry_run" != "true" && -z "${FAVORCLI_FAVOR_CLI_PATH:-}" ]]; then
    echo "Error: favor_cli_path is required when dry_run is false." >&2
    exit 1
  fi

  mkdir -p "$FAVORCLI_OUTPUT_DIR"
}
