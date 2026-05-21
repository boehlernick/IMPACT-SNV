#!/usr/bin/env python3
"""Validation helpers for IMPACT-SNV final GDS files.

Backs the production CLI command:

    impact-snv validate-gds

This module is also imported by impact_snv.gds.finalize after each sample is
finalized. GDS node validation is performed through an R/SeqArray backend
because SeqArray and gdsfmt are the authoritative libraries for reading the GDS
layout produced by this workflow.

Deprecated patho_score/patho_score_calc nodes may be omitted entirely. Their
absence is expected for v1.0.0. If deprecated nodes are present, validation
reports a warning in qc-mode=warn and an error in qc-mode=strict.
"""

from __future__ import annotations

import json
import subprocess
import sys
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Iterable, Optional, Sequence

QC_MODES = {"warn", "strict", "off"}
FINAL_GDS_SUFFIX = "_SNV_IMPACT.gds"
PREPRIORITIZATION_MARKER = ".preprioritization"

REQUIRED_FINAL_NODES = [
    "annotation/info/impact_score",
    "annotation/info/impact_score_calc",
    "annotation/info/tier",
    "annotation/info/scoring_gene",
    "annotation/info/scoring_gene_score",
    "annotation/info/FunctionalAnnotation/VarInfo",
    "annotation/info/FunctionalAnnotation/clnsig",
    "annotation/info/FunctionalAnnotation/bravo_af",
    "annotation/info/FunctionalAnnotation/aloft_prediction",
]

DEPRECATED_NODES = [
    "annotation/info/patho_score",
    "annotation/info/patho_score_calc",
]


@dataclass
class GdsValidationResult:
    """Structured validation result for one final GDS file."""

    gds_path: str
    sample_id: str
    is_valid: bool
    variant_count: int = 0
    sample_count: int = 0
    nonzero_impact_score_variants: Optional[int] = None
    tier_counts: dict[str, int] = field(default_factory=dict)
    missing_nodes: list[str] = field(default_factory=list)
    length_mismatches: dict[str, dict[str, int]] = field(default_factory=dict)
    deprecated_nodes_present: list[str] = field(default_factory=list)
    varinfo_preview: list[str] = field(default_factory=list)
    warnings: list[dict[str, str]] = field(default_factory=list)
    errors: list[dict[str, str]] = field(default_factory=list)

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


def infer_sample_id(gds_path: Path) -> str:
    """Infer sample ID from <sample_id>_SNV_IMPACT.gds."""
    name = gds_path.name
    if name.endswith(FINAL_GDS_SUFFIX):
        return name[: -len(FINAL_GDS_SUFFIX)]
    return gds_path.stem


def discover_final_gds_files(input_dir: Path) -> list[Path]:
    """Find final *_SNV_IMPACT.gds files, excluding preprioritization files."""
    files = sorted(input_dir.glob(f"*{FINAL_GDS_SUFFIX}"))
    return [p for p in files if PREPRIORITIZATION_MARKER not in p.name]


def _r_string(value: str) -> str:
    """Return an R-safe quoted string using JSON escaping rules."""
    return json.dumps(value)


def _r_vector(values: Sequence[str]) -> str:
    """Return an R c(...) character vector expression."""
    return "c(" + ", ".join(json.dumps(v) for v in values) + ")"


def _build_r_validation_code(gds_path: Path) -> str:
    """Build R code that validates one GDS and emits a single JSON object.

    Important: json_escape() uses fixed=TRUE for backslash, quote, newline,
    carriage-return, and tab replacements. This avoids the TRE regex
    "Trailing backslash" failure caused by using a literal backslash as a regex.
    """
    gds_r = _r_string(str(gds_path))
    sample_r = _r_string(infer_sample_id(gds_path))
    required_r = _r_vector(REQUIRED_FINAL_NODES)
    deprecated_r = _r_vector(DEPRECATED_NODES)

    return f"""
suppressPackageStartupMessages({{
  if (!requireNamespace('SeqArray', quietly = TRUE)) stop('R package SeqArray is required')
  if (!requireNamespace('gdsfmt', quietly = TRUE)) stop('R package gdsfmt is required')
  library(SeqArray)
  library(gdsfmt)
}})

json_escape <- function(x) {{
  x <- as.character(x)
  x[is.na(x)] <- ''
  x <- gsub('\\\\', '\\\\\\\\', x, fixed = TRUE)
  x <- gsub('"', '\\\\"', x, fixed = TRUE)
  x <- gsub('\n', '\\\\n', x, fixed = TRUE)
  x <- gsub('\r', '\\\\r', x, fixed = TRUE)
  x <- gsub('\t', '\\\\t', x, fixed = TRUE)
  paste0('"', x, '"')
}}

json_value <- function(x) {{
  if (is.null(x)) return('null')
  if (is.logical(x) && length(x) == 1) return(ifelse(isTRUE(x), 'true', 'false'))
  if (is.numeric(x) && length(x) == 1) {{
    if (is.na(x) || is.nan(x) || is.infinite(x)) return('null')
    return(as.character(x))
  }}
  if (is.character(x) && length(x) == 1) return(json_escape(x))
  if (is.atomic(x)) return(paste0('[', paste(vapply(as.list(x), json_value, character(1)), collapse=','), ']'))
  if (is.list(x)) {{
    nms <- names(x)
    if (is.null(nms)) return(paste0('[', paste(vapply(x, json_value, character(1)), collapse=','), ']'))
    parts <- character(length(x))
    for (i in seq_along(x)) {{
      parts[i] <- paste0(json_escape(nms[i]), ':', json_value(x[[i]]))
    }}
    return(paste0('{{', paste(parts, collapse=','), '}}'))
  }}
  json_escape(as.character(x))
}}

make_issue <- function(code, message, severity='error') {{
  list(code=code, severity=severity, message=message)
}}

result <- list(
  gds_path = {gds_r},
  sample_id = {sample_r},
  is_valid = TRUE,
  variant_count = 0,
  sample_count = 0,
  nonzero_impact_score_variants = NULL,
  tier_counts = list(),
  missing_nodes = character(0),
  length_mismatches = list(),
  deprecated_nodes_present = character(0),
  varinfo_preview = character(0),
  warnings = list(),
  errors = list()
)

required_nodes <- {required_r}
deprecated_nodes <- {deprecated_r}

g <- NULL
tryCatch({{
  g <- SeqArray::seqOpen({gds_r})
  variant_id <- SeqArray::seqGetData(g, 'variant.id')
  n <- length(variant_id)
  result$variant_count <- n
  result$sample_count <- length(SeqArray::seqGetData(g, 'sample.id'))

  if (n <= 0) result$errors[[length(result$errors)+1]] <- make_issue('NO_VARIANTS', 'GDS contains zero variants')
  if (result$sample_count <= 0) result$errors[[length(result$errors)+1]] <- make_issue('NO_SAMPLES', 'GDS contains zero samples')

  for (node in required_nodes) {{
    val <- tryCatch(SeqArray::seqGetData(g, node), error = function(e) NULL)
    if (is.null(val)) {{
      result$missing_nodes <- c(result$missing_nodes, node)
      result$errors[[length(result$errors)+1]] <- make_issue('MISSING_NODE', paste0('Missing required node: ', node))
    }} else if (length(val) != n) {{
      result$length_mismatches[[node]] <- list(observed=length(val), expected=n)
      result$errors[[length(result$errors)+1]] <- make_issue('NODE_LENGTH_MISMATCH', paste0('Node ', node, ' has length ', length(val), ', expected ', n))
    }}
  }}

  for (node in deprecated_nodes) {{
    exists <- !is.null(gdsfmt::index.gdsn(g, node, silent=TRUE))
    if (exists) result$deprecated_nodes_present <- c(result$deprecated_nodes_present, node)
  }}

  impact_score <- tryCatch(SeqArray::seqGetData(g, 'annotation/info/impact_score'), error = function(e) NULL)
  if (!is.null(impact_score)) result$nonzero_impact_score_variants <- sum(impact_score > 0, na.rm=TRUE)

  tier <- tryCatch(SeqArray::seqGetData(g, 'annotation/info/tier'), error = function(e) NULL)
  if (!is.null(tier)) {{
    tt <- table(tier, useNA='ifany')
    result$tier_counts <- as.list(as.integer(tt))
    names(result$tier_counts) <- names(tt)
  }}

  varinfo <- tryCatch(SeqArray::seqGetData(g, 'annotation/info/FunctionalAnnotation/VarInfo'), error = function(e) NULL)
  if (!is.null(varinfo)) {{
    result$varinfo_preview <- as.character(utils::head(varinfo, 10))
    bad_varinfo <- is.na(varinfo) | varinfo == '' | !grepl('^[^[:space:]-]+-[0-9]+-[ACGTNacgtn]+-.+$', varinfo)
    if (any(bad_varinfo)) {{
      idx <- which(bad_varinfo)[1]
      result$warnings[[length(result$warnings)+1]] <- make_issue('VARINFO_FORMAT_WARNING', paste0('VarInfo has non-canonical values; first bad index: ', idx), 'warning')
    }}
  }}
}}, error = function(e) {{
  result$errors[[length(result$errors)+1]] <<- make_issue('GDS_VALIDATION_EXCEPTION', conditionMessage(e))
}}, finally = {{
  if (!is.null(g)) try(SeqArray::seqClose(g), silent=TRUE)
}})

result$is_valid <- length(result$errors) == 0
cat(json_value(result), '\n')
"""


def _run_r_validation(gds_path: Path, rscript: str = "Rscript") -> dict[str, Any]:
    code = _build_r_validation_code(gds_path)
    proc = subprocess.run([rscript, "-e", code], text=True, capture_output=True)
    if proc.returncode != 0:
        raise RuntimeError(
            f"R validation backend failed for {gds_path}\nSTDOUT:\n{proc.stdout}\nSTDERR:\n{proc.stderr}"
        )
    stdout = proc.stdout.strip()
    if not stdout:
        raise RuntimeError(f"R validation backend returned no output for {gds_path}")
    line = [x for x in stdout.splitlines() if x.strip()][-1]
    return json.loads(line)


def validate_gds_file(
    gds_path: Path | str,
    *,
    rscript: str = "Rscript",
    qc_mode: str = "warn",
) -> GdsValidationResult:
    """Validate one final IMPACT-SNV GDS file."""
    path = Path(gds_path)
    if qc_mode not in QC_MODES:
        raise ValueError(f"Invalid qc_mode {qc_mode!r}; expected one of {sorted(QC_MODES)}")

    if not path.exists():
        return GdsValidationResult(
            gds_path=str(path),
            sample_id=infer_sample_id(path),
            is_valid=False,
            errors=[{"code": "FILE_NOT_FOUND", "severity": "error", "message": f"GDS file not found: {path}"}],
        )

    data = _run_r_validation(path, rscript=rscript)
    result = GdsValidationResult(
        gds_path=data.get("gds_path", str(path)),
        sample_id=data.get("sample_id", infer_sample_id(path)),
        is_valid=bool(data.get("is_valid", False)),
        variant_count=int(data.get("variant_count") or 0),
        sample_count=int(data.get("sample_count") or 0),
        nonzero_impact_score_variants=data.get("nonzero_impact_score_variants"),
        tier_counts={str(k): int(v) for k, v in (data.get("tier_counts") or {}).items()},
        missing_nodes=list(data.get("missing_nodes") or []),
        length_mismatches=dict(data.get("length_mismatches") or {}),
        deprecated_nodes_present=list(data.get("deprecated_nodes_present") or []),
        varinfo_preview=list(data.get("varinfo_preview") or []),
        warnings=list(data.get("warnings") or []),
        errors=list(data.get("errors") or []),
    )

    if result.deprecated_nodes_present and qc_mode != "off":
        issue = {
            "code": "DEPRECATED_NODES_PRESENT",
            "severity": "warning" if qc_mode == "warn" else "error",
            "message": "Deprecated patho_* nodes are present: " + ", ".join(result.deprecated_nodes_present),
        }
        if qc_mode == "strict":
            result.errors.append(issue)
            result.is_valid = False
        else:
            result.warnings.append(issue)

    if qc_mode == "off":
        result.is_valid = len(result.errors) == 0

    return result


def validate_gds_files(
    gds_files: Iterable[Path | str],
    *,
    rscript: str = "Rscript",
    qc_mode: str = "warn",
) -> list[GdsValidationResult]:
    return [validate_gds_file(Path(p), rscript=rscript, qc_mode=qc_mode) for p in gds_files]


def validate_gds_dir(
    input_dir: Path | str,
    *,
    rscript: str = "Rscript",
    qc_mode: str = "warn",
) -> list[GdsValidationResult]:
    return validate_gds_files(discover_final_gds_files(Path(input_dir)), rscript=rscript, qc_mode=qc_mode)


def print_validation_summary(results: Sequence[GdsValidationResult]) -> None:
    """Print a concise human-readable validation summary."""
    for res in results:
        status = "OK" if res.is_valid else "FAIL"
        print(f"[{status}] {res.sample_id}: {res.gds_path}")
        print(f"  variants: {res.variant_count}")
        print(f"  samples: {res.sample_count}")
        if res.nonzero_impact_score_variants is not None:
            print(f"  nonzero impact_score variants: {res.nonzero_impact_score_variants}")
        if res.tier_counts:
            print(f"  tier counts: {res.tier_counts}")
        if res.varinfo_preview:
            print(f"  VarInfo preview: {', '.join(res.varinfo_preview[:3])}")
        for warning in res.warnings:
            print(f"  WARNING [{warning.get('code')}]: {warning.get('message')}")
        for error in res.errors:
            print(f"  ERROR [{error.get('code')}]: {error.get('message')}")


def write_validation_manifest(results: Sequence[GdsValidationResult], path: Path | str) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "command": "impact-snv validate-gds",
        "result_count": len(results),
        "valid_count": sum(1 for r in results if r.is_valid),
        "invalid_count": sum(1 for r in results if not r.is_valid),
        "results": [r.to_dict() for r in results],
    }
    path.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")


def run_validate_gds(args: Any) -> int:
    """CLI adapter used by impact_snv.cli."""
    qc_mode = getattr(args, "qc_mode", "warn")
    rscript = getattr(args, "rscript", "Rscript")

    if getattr(args, "gds", None):
        results = [validate_gds_file(args.gds, rscript=rscript, qc_mode=qc_mode)]
    elif getattr(args, "input_dir", None):
        files = discover_final_gds_files(Path(args.input_dir))
        if not files:
            print(f"No final *{FINAL_GDS_SUFFIX} files found in {args.input_dir}", file=sys.stderr)
            return 1
        results = validate_gds_files(files, rscript=rscript, qc_mode=qc_mode)
    else:
        print("Either --gds or --input-dir is required", file=sys.stderr)
        return 2

    print_validation_summary(results)
    return 1 if any(not r.is_valid for r in results) else 0


__all__ = [
    "DEPRECATED_NODES",
    "FINAL_GDS_SUFFIX",
    "GdsValidationResult",
    "REQUIRED_FINAL_NODES",
    "discover_final_gds_files",
    "infer_sample_id",
    "print_validation_summary",
    "run_validate_gds",
    "validate_gds_dir",
    "validate_gds_file",
    "validate_gds_files",
    "write_validation_manifest",
]
