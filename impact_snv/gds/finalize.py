#!/usr/bin/env python3
"""Finalize IMPACT-SNV pre-prioritization GDS files.

This module backs the production CLI command:

    impact-snv finalize-gds

The command converts per-sample pre-prioritization GDS files into final
IMPACT-VIS-ready <sample_id>_SNV_IMPACT.gds files by orchestrating two R
resources packaged with IMPACT-SNV:

  1. impact_prioritize_gds.R
     Adds IMPACT scores, tier labels, scoring gene fields, and ClinVar flags.

  2. add_impact_vis_compat_nodes.R
     Adds IMPACT-VIS compatibility nodes, especially
     annotation/info/FunctionalAnnotation/VarInfo in chr-pos-ref-alt format.

After each final file is produced, this module validates it using
impact_snv.gds.validate.validate_gds_file and writes a finalize_manifest.json.
"""

from __future__ import annotations

import json
import shutil
import subprocess
import sys
from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone
from importlib.resources import files
from pathlib import Path
from typing import Any, Optional, Sequence

from impact_snv.gds.validate import GdsValidationResult, validate_gds_file


PREPRIORITIZATION_SUFFIX = "_SNV_IMPACT.preprioritization.gds"
FINAL_SUFFIX = "_SNV_IMPACT.gds"
VERSION = "1.0.0a0"


@dataclass
class CommandResult:
    """Captured subprocess result for manifest/debugging."""

    command: list[str]
    returncode: int
    stdout_log: str
    stderr_log: str


@dataclass
class FinalizeSampleResult:
    """Structured result for one sample finalization."""

    sample_id: str
    input_gds: str
    scored_intermediate_gds: str
    output_gds: str
    status: str
    message: str = ""
    validation: Optional[dict[str, Any]] = None
    warnings: list[dict[str, str]] = field(default_factory=list)
    errors: list[dict[str, str]] = field(default_factory=list)
    commands: list[CommandResult] = field(default_factory=list)

    def to_dict(self) -> dict[str, Any]:
        d = asdict(self)
        return d


def utc_now_iso() -> str:
    return datetime.now(timezone.utc).isoformat()


def infer_sample_id_from_preprior_gds(path: Path) -> str:
    name = path.name
    if not name.endswith(PREPRIORITIZATION_SUFFIX):
        raise ValueError(f"Not a pre-prioritization GDS file: {path}")
    return name[: -len(PREPRIORITIZATION_SUFFIX)]


def discover_preprioritization_gds_files(input_dir: Path) -> list[Path]:
    """Discover *_SNV_IMPACT.preprioritization.gds files."""
    return sorted(input_dir.glob(f"*{PREPRIORITIZATION_SUFFIX}"))


def resource_script(name: str) -> Path:
    """Resolve an R resource script packaged under impact_snv/resources."""
    candidate = files("impact_snv.resources").joinpath(name)
    path = Path(str(candidate))
    if not path.exists():
        raise FileNotFoundError(
            f"Required resource script not found: {name}. "
            "Ensure pyproject.toml includes package-data for impact_snv/resources/*.R."
        )
    return path


def issue(code: str, message: str, severity: str = "warning") -> dict[str, str]:
    return {"code": code, "severity": severity, "message": message}


def run_command(
    command: Sequence[str],
    *,
    log_prefix: Path,
    label: str,
) -> CommandResult:
    """Run a command, write stdout/stderr logs, and return metadata."""
    log_prefix.parent.mkdir(parents=True, exist_ok=True)
    stdout_log = log_prefix.with_suffix(f".{label}.stdout.log")
    stderr_log = log_prefix.with_suffix(f".{label}.stderr.log")

    proc = subprocess.run(
        list(map(str, command)),
        text=True,
        capture_output=True,
    )
    stdout_log.write_text(proc.stdout or "", encoding="utf-8")
    stderr_log.write_text(proc.stderr or "", encoding="utf-8")

    result = CommandResult(
        command=list(map(str, command)),
        returncode=proc.returncode,
        stdout_log=str(stdout_log),
        stderr_log=str(stderr_log),
    )

    if proc.returncode != 0:
        raise RuntimeError(
            f"Command failed with exit code {proc.returncode}: {' '.join(map(str, command))}\n"
            f"STDOUT log: {stdout_log}\n"
            f"STDERR log: {stderr_log}\n"
            f"STDERR preview:\n{(proc.stderr or '')[-4000:]}"
        )

    return result


def build_prioritize_command(
    *,
    rscript: str,
    prioritize_script: Path,
    input_gds: Path,
    gene_list: Path,
    output_gds: Path,
    force: bool,
) -> list[str]:
    cmd = [
        rscript,
        str(prioritize_script),
        "--input-gds",
        str(input_gds),
        "--gene-list",
        str(gene_list),
        "--output-gds",
        str(output_gds),
    ]
    if force:
        cmd.append("--force")
    return cmd


def build_compat_command(
    *,
    rscript: str,
    compat_script: Path,
    input_gds: Path,
    output_gds: Path,
    force: bool,
    no_optimize: bool,
) -> list[str]:
    cmd = [
        rscript,
        str(compat_script),
        "--input-gds",
        str(input_gds),
        "--output-gds",
        str(output_gds),
    ]
    if force:
        cmd.append("--force")
    if no_optimize:
        cmd.append("--no-optimize")
    return cmd


def finalize_one_sample(
    *,
    sample_id: str,
    input_gds: Path,
    out_dir: Path,
    gene_list: Path,
    rscript: str,
    prioritize_script: Path,
    compat_script: Path,
    qc_mode: str,
    force: bool,
    no_optimize: bool,
    keep_intermediate: bool,
    logs_dir: Path,
    tmp_dir: Path,
) -> FinalizeSampleResult:
    """Finalize one pre-prioritization GDS file."""
    final_gds = out_dir / f"{sample_id}{FINAL_SUFFIX}"
    scored_gds = tmp_dir / f"{sample_id}_SNV_IMPACT.scored.tmp.gds"
    log_prefix = logs_dir / sample_id

    result = FinalizeSampleResult(
        sample_id=sample_id,
        input_gds=str(input_gds),
        scored_intermediate_gds=str(scored_gds),
        output_gds=str(final_gds),
        status="pending",
    )

    try:
        if final_gds.exists() and not force:
            raise FileExistsError(f"Output exists; use --force to overwrite: {final_gds}")

        # Ensure stale temporary scored output does not block the R script.
        if scored_gds.exists():
            scored_gds.unlink()

        if final_gds.exists() and force:
            final_gds.unlink()

        prioritize_cmd = build_prioritize_command(
            rscript=rscript,
            prioritize_script=prioritize_script,
            input_gds=input_gds,
            gene_list=gene_list,
            output_gds=scored_gds,
            force=True,
        )
        result.commands.append(
            run_command(prioritize_cmd, log_prefix=log_prefix, label="prioritize")
        )

        compat_cmd = build_compat_command(
            rscript=rscript,
            compat_script=compat_script,
            input_gds=scored_gds,
            output_gds=final_gds,
            force=True,
            no_optimize=no_optimize,
        )
        result.commands.append(
            run_command(compat_cmd, log_prefix=log_prefix, label="viscompat")
        )

        validation: GdsValidationResult = validate_gds_file(
            final_gds,
            rscript=rscript,
            qc_mode=qc_mode,
        )
        result.validation = validation.to_dict()
        result.warnings.extend(validation.warnings)
        result.errors.extend(validation.errors)

        if not validation.is_valid:
            result.status = "failed"
            result.message = "Final GDS validation failed"
        else:
            result.status = "ok"
            result.message = ""

    except Exception as exc:
        result.status = "failed"
        result.message = str(exc)
        result.errors.append(issue("FINALIZE_EXCEPTION", str(exc), severity="error"))

    finally:
        if scored_gds.exists() and not keep_intermediate:
            try:
                scored_gds.unlink()
            except Exception as exc:  # pragma: no cover - cleanup warning only
                result.warnings.append(
                    issue(
                        "INTERMEDIATE_CLEANUP_WARNING",
                        f"Could not remove intermediate GDS {scored_gds}: {exc}",
                    )
                )

    return result


def write_manifest(
    *,
    manifest_path: Path,
    input_dir: Path,
    out_dir: Path,
    gene_list: Path,
    qc_mode: str,
    prioritize_script: Path,
    compat_script: Path,
    results: Sequence[FinalizeSampleResult],
) -> None:
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "command": "impact-snv finalize-gds",
        "version": VERSION,
        "created_utc": utc_now_iso(),
        "input_dir": str(input_dir),
        "out_dir": str(out_dir),
        "gene_list": str(gene_list),
        "qc_mode": qc_mode,
        "prioritize_script": str(prioritize_script),
        "compat_script": str(compat_script),
        "result_count": len(results),
        "ok_count": sum(1 for r in results if r.status == "ok"),
        "failed_count": sum(1 for r in results if r.status != "ok"),
        "results": [r.to_dict() for r in results],
    }
    manifest_path.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")


def run_finalize_gds(args: Any) -> int:
    """CLI adapter used by impact_snv.cli for finalize-gds."""
    input_dir = Path(args.input_dir)
    gene_list = Path(args.gene_list)
    out_dir = Path(args.out_dir)
    qc_mode = getattr(args, "qc_mode", "warn")
    samples_filter = set(getattr(args, "samples", None) or [])
    rscript = getattr(args, "rscript", "Rscript")
    force = bool(getattr(args, "force", False))
    no_optimize = bool(getattr(args, "no_optimize", False))
    keep_intermediate = bool(getattr(args, "keep_intermediate", False))
    manifest_json = getattr(args, "manifest_json", None)

    out_dir.mkdir(parents=True, exist_ok=True)
    logs_dir = out_dir / "logs"
    tmp_dir = out_dir / ".tmp"
    logs_dir.mkdir(parents=True, exist_ok=True)
    tmp_dir.mkdir(parents=True, exist_ok=True)

    prioritize_script = resource_script("impact_prioritize_gds.R")
    compat_script = resource_script("add_impact_vis_compat_nodes.R")

    preprior_files = discover_preprioritization_gds_files(input_dir)
    if samples_filter:
        preprior_files = [
            p for p in preprior_files if infer_sample_id_from_preprior_gds(p) in samples_filter
        ]

    if not preprior_files:
        print(
            f"No pre-prioritization GDS files found in {input_dir}"
            + (f" for samples: {', '.join(sorted(samples_filter))}" if samples_filter else ""),
            file=sys.stderr,
        )
        return 1

    results: list[FinalizeSampleResult] = []
    for input_gds in preprior_files:
        sample_id = infer_sample_id_from_preprior_gds(input_gds)
        print(f"Finalizing {sample_id}: {input_gds}", flush=True)
        sample_result = finalize_one_sample(
            sample_id=sample_id,
            input_gds=input_gds,
            out_dir=out_dir,
            gene_list=gene_list,
            rscript=rscript,
            prioritize_script=prioritize_script,
            compat_script=compat_script,
            qc_mode=qc_mode,
            force=force,
            no_optimize=no_optimize,
            keep_intermediate=keep_intermediate,
            logs_dir=logs_dir,
            tmp_dir=tmp_dir,
        )
        results.append(sample_result)

        if sample_result.status == "ok":
            validation = sample_result.validation or {}
            print(
                f"  OK: {sample_result.output_gds} "
                f"variants={validation.get('variant_count', 'NA')} "
                f"nonzero={validation.get('nonzero_impact_score_variants', 'NA')}",
                flush=True,
            )
        else:
            print(f"  FAILED: {sample_result.message}", file=sys.stderr, flush=True)

    manifest_path = Path(manifest_json) if manifest_json else out_dir / "finalize_manifest.json"
    write_manifest(
        manifest_path=manifest_path,
        input_dir=input_dir,
        out_dir=out_dir,
        gene_list=gene_list,
        qc_mode=qc_mode,
        prioritize_script=prioritize_script,
        compat_script=compat_script,
        results=results,
    )
    print(f"Manifest written: {manifest_path}", flush=True)

    failed = [r for r in results if r.status != "ok"]
    if failed:
        print(f"Finalization completed with {len(failed)} failed sample(s).", file=sys.stderr)
        return 1

    print("Finalization completed successfully.", flush=True)

    # Remove empty tmp dir unless user kept intermediates.
    if not keep_intermediate:
        try:
            if tmp_dir.exists() and not any(tmp_dir.iterdir()):
                tmp_dir.rmdir()
        except Exception:
            pass

    return 0


__all__ = [
    "FINAL_SUFFIX",
    "PREPRIORITIZATION_SUFFIX",
    "FinalizeSampleResult",
    "discover_preprioritization_gds_files",
    "infer_sample_id_from_preprior_gds",
    "run_finalize_gds",
]
