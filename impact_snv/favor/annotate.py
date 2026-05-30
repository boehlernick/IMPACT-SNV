#!/usr/bin/env python3
"""Step 3 annotation backend wrappers for IMPACT-SNV.

`impact-snv favor-annotate` now supports a backend selector:
- `favor-cli` (default): invokes `favor annotate`.
- `legacy-favorannotator`: compatibility wrapper that stages existing legacy
  FAVORannotator-style outputs into the canonical `<out_prefix>.annotated` /
  `<out_prefix>.genotypes` layout without changing Step 4 contracts.
"""

from __future__ import annotations

import shutil
import subprocess
import os
from pathlib import Path
from typing import Any, Optional

from impact_snv.favor.common import (
    FavorRunResult,
    chromosome_partition_dirs,
    ensure_out_dir,
    observed_mtimes,
    observed_sizes,
    remove_existing_outputs,
    resolve_favor_bin,
    run_logged_command_with_progress,
    utc_now_iso,
    validate_reference_build,
    write_manifest,
)

BACKEND_FAVOR_CLI = "favor-cli"
BACKEND_LEGACY = "legacy-favorannotator"
BACKEND_FAVOR_CLI_SKELETON = "favor-cli-skeleton"
SUPPORTED_BACKENDS = (BACKEND_FAVOR_CLI, BACKEND_LEGACY, BACKEND_FAVOR_CLI_SKELETON)


def infer_prefix_from_ingested_dir(path: Path) -> str:
    return path.name[: -len(".ingested")] if path.name.endswith(".ingested") else path.name


def _warning(code: str, message: str) -> dict[str, str]:
    return {"code": code, "severity": "warning", "message": message}


def _required_dir(value: Optional[str | Path], field_name: str) -> Path:
    if value is None:
        raise ValueError(f"{field_name} is required for selected backend")
    path = Path(value).expanduser().resolve()
    if not path.exists() or not path.is_dir():
        raise FileNotFoundError(f"Directory not found for {field_name}: {path}")
    return path


def _required_file(value: Optional[str | Path], field_name: str) -> Path:
    if value is None:
        raise ValueError(f"{field_name} is required for selected backend")
    path = Path(value).expanduser().resolve()
    if not path.exists() or not path.is_file():
        raise FileNotFoundError(f"File not found for {field_name}: {path}")
    return path


def _remove_path(path: Path) -> None:
    if path.is_symlink() or path.is_file():
        path.unlink()
    elif path.is_dir():
        shutil.rmtree(path)


def _stage_legacy_output(src: Path, dst: Path, *, force: bool, mode: str) -> str:
    """Stage a legacy directory into canonical output layout."""
    src = src.resolve()
    dst = dst.resolve()
    if src == dst:
        return "in-place"

    if dst.exists():
        if dst.resolve() == src:
            return "already-staged"
        if not force:
            raise FileExistsError(f"Destination exists; use --force to overwrite: {dst}")
        _remove_path(dst)

    dst.parent.mkdir(parents=True, exist_ok=True)

    if mode == "none":
        return "none"
    if mode == "copy":
        shutil.copytree(src, dst)
        return "copy"

    # default: symlink; fall back to copy if symlink fails in the environment
    try:
        dst.symlink_to(src, target_is_directory=True)
        return "symlink"
    except OSError:
        shutil.copytree(src, dst)
        return "copy-fallback"


def run_favor_cli_backend(args: Any) -> int:
    ingested_dir = _required_dir(getattr(args, "ingested_dir", None), "--ingested-dir")
    reference_build = getattr(args, "reference_build", "GRCh38")
    validate_reference_build(reference_build)

    favor_bin = resolve_favor_bin(getattr(args, "favor_bin", "favor"))
    out_dir = ensure_out_dir(Path(args.out_dir))
    out_prefix = getattr(args, "out_prefix", None) or infer_prefix_from_ingested_dir(ingested_dir)
    force = bool(getattr(args, "force", False))

    progress_interval = int(getattr(args, "progress_interval_seconds", 60) or 60)
    quiet = bool(getattr(args, "quiet", False))
    tail_log_lines = int(getattr(args, "tail_log_lines", 0) or 0)
    max_chars = int(getattr(args, "max_progress_log_line_chars", 300) or 300)
    progress_mode = getattr(args, "progress_mode", "normal") or "normal"

    annotated_dir = out_dir / f"{out_prefix}.annotated"
    genotypes_dir = out_dir / f"{out_prefix}.genotypes"
    samples_file = genotypes_dir / "samples.txt"
    removed_paths = remove_existing_outputs([annotated_dir], force=force)

    try:
        ingested_arg = str(ingested_dir.relative_to(out_dir))
    except ValueError:
        ingested_arg = str(ingested_dir)

    stdout_log = out_dir / "logs" / f"{out_prefix}.favor_annotate.stdout.log"
    stderr_log = out_dir / "logs" / f"{out_prefix}.favor_annotate.stderr.log"
    command = [favor_bin, "annotate", ingested_arg]

    start = utc_now_iso()
    rc = run_logged_command_with_progress(
        command,
        cwd=out_dir,
        stdout_log=stdout_log,
        stderr_log=stderr_log,
        progress_label="favor-annotate",
        progress_interval_seconds=progress_interval,
        quiet=quiet,
        watch_paths=[annotated_dir],
        tail_log_lines=tail_log_lines,
        max_progress_log_line_chars=max_chars,
        progress_mode=progress_mode,
    )
    end = utc_now_iso()

    ann_parts = chromosome_partition_dirs(annotated_dir)
    gt_parts = chromosome_partition_dirs(genotypes_dir)
    status = "ok"
    messages: list[str] = []
    warnings: list[dict[str, str]] = []

    if rc != 0:
        status = "failed"
        messages.append(f"favor annotate exited with return code {rc}; see {stderr_log}")
    if not annotated_dir.exists():
        status = "failed"
        messages.append(f"Expected annotated output was not created: {annotated_dir}")
    if not ann_parts:
        status = "failed"
        messages.append("No chromosome partition directories were found in annotated output")
    if not genotypes_dir.exists():
        warnings.append(
            _warning(
                "GENOTYPES_OUTPUT_NOT_CREATED",
                f"favor annotate did not create {genotypes_dir}; genotype extraction is handled by impact-snv extract-genotypes.",
            )
        )
    if not samples_file.exists():
        warnings.append(
            _warning(
                "GENOTYPE_SAMPLES_FILE_NOT_CREATED",
                f"favor annotate did not create {samples_file}; genotype extraction is handled by impact-snv extract-genotypes.",
            )
        )

    msg = "; ".join(messages)
    result = FavorRunResult(
        command="favor-annotate",
        favor_command=list(map(str, command)),
        cwd=str(out_dir),
        out_dir=str(out_dir),
        out_prefix=out_prefix,
        reference_build=reference_build,
        favor_bin=getattr(args, "favor_bin", "favor"),
        force=force,
        return_code=rc,
        start_time_utc=start,
        end_time_utc=end,
        stdout_log=str(stdout_log),
        stderr_log=str(stderr_log),
        progress_interval_seconds=progress_interval,
        quiet=quiet,
        tail_log_lines=tail_log_lines,
        max_progress_log_line_chars=max_chars,
        progress_mode=progress_mode,
        removed_paths=removed_paths,
        observed_output_sizes_bytes=observed_sizes({"annotated_dir": annotated_dir}),
        observed_output_latest_mtime_utc=observed_mtimes({"annotated_dir": annotated_dir}),
        status=status,
        message=msg,
    )
    payload = result.to_dict()
    payload.update(
        {
            "backend": BACKEND_FAVOR_CLI,
            "ingested_dir": str(ingested_dir),
            "annotated_dir": str(annotated_dir.resolve()),
            "genotypes_dir": str(genotypes_dir.resolve()),
            "samples_file": str(samples_file.resolve()),
            "annotated_chromosome_partitions": ann_parts,
            "genotype_chromosome_partitions": gt_parts,
            "warnings": warnings,
        }
    )
    manifest = Path(getattr(args, "manifest_json", None) or (out_dir / f"{out_prefix}.annotate_manifest.json"))
    write_manifest(manifest, payload)

    if status != "ok":
        print(msg)
        return 1

    if warnings:
        for warning in warnings:
            print(f"WARNING [{warning['code']}]: {warning['message']}")
    print(f"FAVOR annotate completed successfully: {annotated_dir}")
    print(f"Manifest written: {manifest}")
    return 0


def run_legacy_favorannotator_backend(args: Any) -> int:
    """Compatibility backend for legacy FAVORannotator outputs.

    This backend does not invoke legacy DNAnexus applets directly. It wraps
    existing legacy-generated directories and stages them into canonical output
    layout for downstream build/finalize steps.
    """

    out_dir = ensure_out_dir(Path(args.out_dir))
    reference_build = getattr(args, "reference_build", "GRCh38")
    validate_reference_build(reference_build)

    source_annotated = _required_dir(
        getattr(args, "legacy_annotated_dir", None),
        "--legacy-annotated-dir",
    )
    source_genotypes = None
    legacy_genotypes = getattr(args, "legacy_genotypes_dir", None)
    if legacy_genotypes is not None:
        source_genotypes = _required_dir(legacy_genotypes, "--legacy-genotypes-dir")

    out_prefix = getattr(args, "out_prefix", None) or source_annotated.name.replace(".annotated", "")
    force = bool(getattr(args, "force", False))
    stage_mode = getattr(args, "legacy_stage_mode", "symlink")
    if stage_mode not in {"symlink", "copy", "none"}:
        raise ValueError(f"Unsupported --legacy-stage-mode {stage_mode!r}")

    target_annotated = out_dir / f"{out_prefix}.annotated"
    target_genotypes = out_dir / f"{out_prefix}.genotypes"

    stage_ann = _stage_legacy_output(source_annotated, target_annotated, force=force, mode=stage_mode)
    stage_gt = "none"
    if source_genotypes is not None:
        stage_gt = _stage_legacy_output(source_genotypes, target_genotypes, force=force, mode=stage_mode)

    ann_parts = chromosome_partition_dirs(target_annotated if stage_mode != "none" else source_annotated)
    gt_parts = chromosome_partition_dirs(target_genotypes) if source_genotypes is not None and stage_mode != "none" else []

    warnings: list[dict[str, str]] = []
    if not ann_parts:
        warnings.append(
            _warning(
                "LEGACY_ANNOTATED_PARTITIONS_MISSING",
                f"No chromosome= partitions found under {target_annotated if stage_mode != 'none' else source_annotated}.",
            )
        )
    if source_genotypes is None:
        warnings.append(
            _warning(
                "LEGACY_GENOTYPES_NOT_PROVIDED",
                "--legacy-genotypes-dir was not provided; run impact-snv extract-genotypes separately for downstream build-gds.",
            )
        )

    manifest = Path(getattr(args, "manifest_json", None) or (out_dir / f"{out_prefix}.annotate_manifest.json"))
    payload = {
        "command": "impact-snv favor-annotate",
        "backend": BACKEND_LEGACY,
        "status": "ok",
        "message": "Legacy FAVORannotator compatibility backend staged existing outputs.",
        "reference_build": reference_build,
        "out_dir": str(out_dir),
        "out_prefix": out_prefix,
        "legacy_annotated_dir": str(source_annotated),
        "legacy_genotypes_dir": str(source_genotypes) if source_genotypes is not None else None,
        "annotated_dir": str(target_annotated),
        "genotypes_dir": str(target_genotypes),
        "annotated_stage_mode": stage_ann,
        "genotypes_stage_mode": stage_gt,
        "annotated_chromosome_partitions": ann_parts,
        "genotype_chromosome_partitions": gt_parts,
        "warnings": warnings,
        "created_utc": utc_now_iso(),
        "observed_output_sizes_bytes": observed_sizes({
            "annotated_dir": target_annotated if stage_mode != "none" else source_annotated,
            "genotypes_dir": target_genotypes if source_genotypes is not None and stage_mode != "none" else (source_genotypes or out_dir),
        }),
        "observed_output_latest_mtime_utc": observed_mtimes({
            "annotated_dir": target_annotated if stage_mode != "none" else source_annotated,
            "genotypes_dir": target_genotypes if source_genotypes is not None and stage_mode != "none" else (source_genotypes or out_dir),
        }),
    }
    write_manifest(manifest, payload)

    for warning in warnings:
        print(f"WARNING [{warning['code']}]: {warning['message']}")
    print("Legacy FAVORannotator compatibility backend completed.")
    print(f"Manifest written: {manifest}")
    return 0


def run_favor_cli_skeleton_backend(args: Any) -> int:
    """Run step3_favorcli_annotation scaffold in dry-run mode.

    This backend intentionally supports dry-run only and does not mutate GDS or
    produce Step 4-ready annotated outputs.
    """

    reference_build = getattr(args, "reference_build", "GRCh38")
    validate_reference_build(reference_build)
    out_dir = ensure_out_dir(Path(args.out_dir))

    input_file = _required_file(getattr(args, "skeleton_input_file", None), "--skeleton-input-file")
    input_type = getattr(args, "skeleton_input_type", "gds")
    if input_type not in {"gds", "vcf"}:
        raise ValueError(f"Unsupported --skeleton-input-type {input_type!r}")

    # Dry-run only for now to preserve scaffold behavior and avoid guessing CLI schema.
    if not bool(getattr(args, "skeleton_dry_run", True)):
        raise ValueError("favor-cli-skeleton backend currently supports dry-run only")

    out_prefix = getattr(args, "out_prefix", None) or input_file.stem
    skeleton_dir = out_dir / f"{out_prefix}.favorcli_skeleton"
    skeleton_dir.mkdir(parents=True, exist_ok=True)

    repo_root = Path(__file__).resolve().parents[2]
    script_dir = repo_root / "step3_favorcli_annotation" / "src"
    run_script = script_dir / "run_favorcli.sh"
    if not run_script.exists():
        raise FileNotFoundError(f"Missing scaffold script: {run_script}")

    env = os.environ.copy()
    env["FAVORCLI_GDS_FILE"] = str(input_file) if input_type == "gds" else ""
    env["FAVORCLI_VCF_FILE"] = str(input_file) if input_type == "vcf" else ""
    env["FAVORCLI_DATABASE_FILE"] = str(getattr(args, "favor_database_file", "") or "")
    env["FAVORCLI_OUTPUT_DIR"] = str(skeleton_dir)
    env["FAVORCLI_FAVOR_CLI_PATH"] = str(getattr(args, "favor_bin", "favor"))
    env["FAVORCLI_DATABASE_VERSION"] = str(getattr(args, "favor_database_version", "") or "")
    env["FAVORCLI_REFERENCE_BUILD"] = reference_build
    env["FAVORCLI_THREADS"] = str(getattr(args, "skeleton_threads", 4) or 4)
    env["FAVORCLI_MEMORY_BUDGET_GB"] = str(getattr(args, "skeleton_memory_budget_gb", 8) or 8)
    env["FAVORCLI_DRY_RUN"] = "true"

    subprocess.run(["bash", str(run_script)], cwd=str(script_dir), env=env, check=True)

    manifest = Path(getattr(args, "manifest_json", None) or (out_dir / f"{out_prefix}.annotate_manifest.json"))
    payload = {
        "command": "impact-snv favor-annotate",
        "backend": BACKEND_FAVOR_CLI_SKELETON,
        "status": "ok",
        "message": "favor-cli scaffold dry-run completed",
        "reference_build": reference_build,
        "out_dir": str(out_dir),
        "out_prefix": out_prefix,
        "skeleton_input_file": str(input_file),
        "skeleton_input_type": input_type,
        "skeleton_output_dir": str(skeleton_dir),
        "variant_identity_tsv": str(skeleton_dir / "variant_identity.tsv"),
        "validation_report": str(skeleton_dir / "favorcli_validation_report.txt"),
        "skeleton_run_manifest": str(skeleton_dir / "favorcli_run_manifest.json"),
        "created_utc": utc_now_iso(),
        "warnings": [
            _warning(
                "SCAFFOLD_DRY_RUN_ONLY",
                "favor-cli-skeleton backend is dry-run only and does not produce Step 4-ready annotated outputs.",
            )
        ],
    }
    write_manifest(manifest, payload)
    print("FAVOR-CLI skeleton dry-run completed.")
    print(f"Manifest written: {manifest}")
    return 0


def run_favor_annotate(args: Any) -> int:
    backend = getattr(args, "backend", BACKEND_FAVOR_CLI) or BACKEND_FAVOR_CLI
    if backend not in SUPPORTED_BACKENDS:
        raise ValueError(f"Unsupported backend {backend!r}; expected one of {SUPPORTED_BACKENDS}")
    if backend == BACKEND_LEGACY:
        return run_legacy_favorannotator_backend(args)
    if backend == BACKEND_FAVOR_CLI_SKELETON:
        return run_favor_cli_skeleton_backend(args)
    return run_favor_cli_backend(args)
