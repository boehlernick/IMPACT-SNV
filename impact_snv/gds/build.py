#!/usr/bin/env python3
"""Build pre-prioritization IMPACT-SNV GDS files from FAVOR outputs.

Backs the production CLI command:

    impact-snv build-gds

This command coordinates:
  1. impact_snv/gds/flatten.py for each sample/chromosome
  2. impact_snv/resources/favor_flat_to_seqarray_gds.R for GDS writing
  3. per-sample concatenation of chromosome flat parquet files
  4. build_manifest.json and QC outputs

Output GDS files are pre-prioritization files named:

    <sample_id>_SNV_IMPACT.preprioritization.gds

These are consumed by:

    impact-snv finalize-gds
"""

from __future__ import annotations

import json
import re
import subprocess
import sys
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from importlib.resources import files
from pathlib import Path
from typing import Any, Optional, Sequence

import pandas as pd
import pyarrow.parquet as pq

from impact_snv import __version__ as VERSION
from impact_snv.qc.build_qc import run_build_qc

PREPRIORITIZATION_SUFFIX = "_SNV_IMPACT.preprioritization.gds"


@dataclass
class ChromosomeJobResult:
    sample_id: str
    chromosome: str
    flat_parquet: str
    summary_json: str
    gds_path: Optional[str]
    flat_rows: int
    status: str
    message: str = ""


@dataclass
class SampleMergeResult:
    sample_id: str
    merged_gds: Optional[str]
    merged_flat_parquet: Optional[str]
    input_flat_parquets: list[str]
    input_gds_files: list[str]
    variant_count: int
    status: str
    message: str = ""


def utc_now_iso() -> str:
    return datetime.now(timezone.utc).isoformat()


def parse_chromosomes(values: Sequence[str]) -> list[str]:
    out: list[str] = []
    for value in values:
        value = str(value).strip()
        m = re.fullmatch(r"(\d+)-(\d+)", value)
        if m:
            start, end = int(m.group(1)), int(m.group(2))
            step = 1 if start <= end else -1
            out.extend(str(i) for i in range(start, end + step, step))
        else:
            out.append(re.sub(r"^chr", "", value, flags=re.IGNORECASE))
    seen: set[str] = set()
    final: list[str] = []
    for chrom in out:
        if chrom not in seen:
            seen.add(chrom)
            final.append(chrom)
    return final


def chrom_sort_key(chrom: str) -> tuple[int, Any]:
    c = re.sub(r"^chr", "", str(chrom), flags=re.IGNORECASE)
    if c.isdigit():
        return (0, int(c))
    if c.upper() == "X":
        return (1, 23)
    if c.upper() == "Y":
        return (1, 24)
    if c.upper() in {"M", "MT"}:
        return (1, 25)
    return (2, c)


def load_samples(genotypes_dir: Path, requested: Optional[Sequence[str]], all_samples: bool) -> list[str]:
    samples_path = genotypes_dir / "samples.txt"
    if not samples_path.exists():
        raise FileNotFoundError(f"Missing samples.txt: {samples_path}")
    available = [line.strip() for line in samples_path.read_text(encoding="utf-8").splitlines() if line.strip()]
    if not available:
        raise ValueError(f"No samples found in {samples_path}")
    if all_samples or not requested:
        return available
    missing = sorted(set(requested) - set(available))
    if missing:
        raise ValueError(f"Requested samples not present in samples.txt: {missing}; available={available}")
    return list(requested)


def parquet_num_rows(path: Path) -> int:
    if not path.exists():
        return 0
    return int(pq.ParquetFile(path).metadata.num_rows)


def package_path(module: str) -> Path:
    mod = __import__(module, fromlist=["__file__"])
    return Path(mod.__file__).resolve()


def resource_script(name: str) -> Path:
    candidate = files("impact_snv.resources").joinpath(name)
    path = Path(str(candidate))
    if not path.exists():
        raise FileNotFoundError(
            f"Required resource script not found: {name}. "
            "Ensure pyproject.toml includes package-data for impact_snv/resources/*.R."
        )
    return path


def run_cmd(cmd: Sequence[str], *, dry_run: bool = False) -> None:
    print("+", " ".join(map(str, cmd)), flush=True)
    if not dry_run:
        subprocess.run(list(map(str, cmd)), check=True)


def flatten_one(args: Any, sample_id: str, chrom: str, flat_parquet: Path, preview_csv: Optional[Path], summary_json: Path) -> None:
    if flat_parquet.exists() and not args.force:
        print(f"Skipping existing flat parquet: {flat_parquet}", flush=True)
        return
    flat_parquet.parent.mkdir(parents=True, exist_ok=True)
    summary_json.parent.mkdir(parents=True, exist_ok=True)
    cmd = [
        args.python,
        str(args.flatten_script),
        "--annotated-dir", str(args.annotated_dir),
        "--genotypes-dir", str(args.genotypes_dir),
        "--gene-list", str(args.gene_list),
        "--sample-id", sample_id,
        "--chromosome", chrom,
        "--dosage-threshold", str(args.dosage_threshold),
        "--out", str(flat_parquet),
        "--summary-json", str(summary_json),
    ]
    if preview_csv is not None:
        preview_csv.parent.mkdir(parents=True, exist_ok=True)
        cmd.extend(["--preview-csv", str(preview_csv)])
    run_cmd(cmd, dry_run=args.dry_run)


def write_gds_one(args: Any, sample_id: str, flat_parquet: Path, gds_path: Path) -> None:
    if gds_path.exists() and not args.force:
        print(f"Skipping existing GDS: {gds_path}", flush=True)
        return
    gds_path.parent.mkdir(parents=True, exist_ok=True)
    cmd = [
        args.rscript,
        str(args.gds_writer),
        "--input", str(flat_parquet),
        "--output", str(gds_path),
        "--sample-id", sample_id,
    ]
    if args.force:
        cmd.append("--force")
    if args.no_optimize:
        cmd.append("--no-optimize")
    if args.keep_temp_vcf:
        cmd.extend(["--keep-temp-vcf", "--temp-vcf", str(gds_path.with_suffix(".temp.vcf"))])
    run_cmd(cmd, dry_run=args.dry_run)


def concat_flat_parquets(flat_paths: Sequence[Path], output_path: Path, sample_id: str) -> int:
    frames: list[pd.DataFrame] = []
    for path in flat_paths:
        if path.exists() and parquet_num_rows(path) > 0:
            frames.append(pd.read_parquet(path, engine="pyarrow"))
    if not frames:
        return 0
    df = pd.concat(frames, ignore_index=True)
    df = df[df["sample_id"].astype(str) == sample_id].copy()
    if df.empty:
        return 0
    df["_chrom_rank"] = df["chromosome"].map(lambda x: chrom_sort_key(str(x))[0])
    df["_chrom_value"] = df["chromosome"].map(lambda x: chrom_sort_key(str(x))[1])
    sort_cols = ["_chrom_rank", "_chrom_value", "position", "ref", "alt"]
    df = df.sort_values(sort_cols, kind="mergesort")
    df = df.drop(columns=["_chrom_rank", "_chrom_value"])
    df["variant_id"] = range(1, len(df) + 1)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_parquet(output_path, index=False, engine="pyarrow")
    return int(len(df))


def merge_sample_by_concat(args: Any, sample_id: str, flat_paths: Sequence[Path], chrom_gds_files: Sequence[Path], output_flat: Path, output_gds: Path) -> SampleMergeResult:
    if output_gds.exists() and not args.force:
        return SampleMergeResult(
            sample_id, str(output_gds), str(output_flat), [str(x) for x in flat_paths],
            [str(x) for x in chrom_gds_files], parquet_num_rows(output_flat) if output_flat.exists() else 0,
            "skipped", "Merged GDS already exists"
        )
    try:
        n = concat_flat_parquets(flat_paths, output_flat, sample_id)
        if n == 0:
            return SampleMergeResult(sample_id, None, str(output_flat), [str(x) for x in flat_paths], [str(x) for x in chrom_gds_files], 0, "skipped", "No rows to merge")
        write_gds_one(args, sample_id, output_flat, output_gds)
        return SampleMergeResult(sample_id, str(output_gds), str(output_flat), [str(x) for x in flat_paths], [str(x) for x in chrom_gds_files], n, "ok")
    except Exception as exc:
        return SampleMergeResult(sample_id, None, str(output_flat), [str(x) for x in flat_paths], [str(x) for x in chrom_gds_files], 0, "error", str(exc))


def build_manifest_payload(args: Any, samples: Sequence[str], chromosomes: Sequence[str], jobs: Sequence[ChromosomeJobResult], merges: Sequence[SampleMergeResult]) -> dict[str, Any]:
    return {
        "command": "impact-snv build-gds",
        "version": VERSION,
        "created_utc": utc_now_iso(),
        "annotated_dir": str(args.annotated_dir),
        "genotypes_dir": str(args.genotypes_dir),
        "gene_list": str(args.gene_list),
        "out_dir": str(args.out_dir),
        "samples": list(samples),
        "chromosomes": list(chromosomes),
        "dosage_threshold": args.dosage_threshold,
        "merge_per_sample": args.merge_per_sample,
        "skip_chrom_gds": args.skip_chrom_gds,
        "qc_mode": args.qc_mode,
        "jobs": [asdict(x) for x in jobs],
        "sample_merges": [asdict(x) for x in merges],
    }


def run_build_gds(args: Any) -> int:
    """CLI adapter used by impact_snv.cli for build-gds."""
    args.out_dir.mkdir(parents=True, exist_ok=True)
    flat_dir = args.out_dir / "flat"
    summary_dir = args.out_dir / "summary"
    chrom_gds_dir = args.out_dir / "gds_by_chrom"
    merged_flat_dir = args.out_dir / "flat_merged"
    merged_gds_dir = args.out_dir / "gds_merged"
    for d in [flat_dir, summary_dir, chrom_gds_dir, merged_flat_dir, merged_gds_dir]:
        d.mkdir(parents=True, exist_ok=True)

    samples = load_samples(args.genotypes_dir, args.samples, args.all_samples)
    chromosomes = parse_chromosomes(args.chromosomes)

    jobs: list[ChromosomeJobResult] = []
    sample_flat_paths: dict[str, list[Path]] = {s: [] for s in samples}
    sample_chrom_gds: dict[str, list[Path]] = {s: [] for s in samples}

    for sample_id in samples:
        for chrom in chromosomes:
            flat_parquet = flat_dir / f"{sample_id}_chr{chrom}.flat.parquet"
            preview_csv = flat_dir / f"{sample_id}_chr{chrom}.preview.head.csv" if args.write_preview_csv else None
            summary_json = summary_dir / f"{sample_id}_chr{chrom}.summary.json"
            chrom_gds = chrom_gds_dir / f"{sample_id}_chr{chrom}.gds"
            try:
                flatten_one(args, sample_id, chrom, flat_parquet, preview_csv, summary_json)
                rows = 0 if args.dry_run else parquet_num_rows(flat_parquet)
                sample_flat_paths[sample_id].append(flat_parquet)
                gds_path: Optional[str] = None
                status = "ok"
                message = ""
                if rows == 0:
                    status = "skipped"
                    message = "No rows after filters"
                elif not args.skip_chrom_gds:
                    write_gds_one(args, sample_id, flat_parquet, chrom_gds)
                    gds_path = str(chrom_gds)
                    sample_chrom_gds[sample_id].append(chrom_gds)
                jobs.append(ChromosomeJobResult(sample_id, chrom, str(flat_parquet), str(summary_json), gds_path, rows, status, message))
            except Exception as exc:
                jobs.append(ChromosomeJobResult(sample_id, chrom, str(flat_parquet), str(summary_json), None, 0, "error", str(exc)))
                if args.qc_mode == "strict":
                    raise

    merges: list[SampleMergeResult] = []
    if args.merge_per_sample:
        for sample_id in samples:
            output_flat = merged_flat_dir / f"{sample_id}_all_chromosomes.flat.parquet"
            output_gds = merged_gds_dir / f"{sample_id}{PREPRIORITIZATION_SUFFIX}"
            merges.append(
                merge_sample_by_concat(
                    args, sample_id, sample_flat_paths[sample_id], sample_chrom_gds[sample_id], output_flat, output_gds
                )
            )

    manifest_path = args.manifest_json or (args.out_dir / "build_manifest.json")
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    manifest = build_manifest_payload(args, samples, chromosomes, jobs, merges)
    manifest_path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(f"Build manifest written: {manifest_path}")

    qc_result = run_build_qc(manifest_path, args.out_dir / "qc", qc_mode=args.qc_mode, print_summary=True)
    if not qc_result.is_valid:
        print("Build QC failed in strict mode", file=sys.stderr)
        return 2

    failed_jobs = [j for j in jobs if j.status == "error"]
    failed_merges = [m for m in merges if m.status == "error"]
    if failed_jobs or failed_merges:
        print(f"Build completed with errors: jobs={len(failed_jobs)} merges={len(failed_merges)}", file=sys.stderr)
        return 1

    print("Build completed successfully.")
    return 0


__all__ = [
    "ChromosomeJobResult",
    "SampleMergeResult",
    "parse_chromosomes",
    "run_build_gds",
]
