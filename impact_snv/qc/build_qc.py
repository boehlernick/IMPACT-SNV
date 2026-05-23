#!/usr/bin/env python3
"""Build-manifest QC for IMPACT-SNV GDS construction.

This module backs the planned production CLI command:

    impact-snv qc-build

and is also imported by impact_snv.gds.build after build-gds completes.

Default v1.0.0 policy:
  - zero-row autosomes are warnings, not fatal failures
  - qc-mode=strict may convert zero-row autosomes into validation errors
"""

from __future__ import annotations

import argparse
import json
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Any, Optional

import pandas as pd

QC_MODES = {"warn", "strict", "off"}


@dataclass
class BuildQcResult:
    manifest: str
    out_dir: str
    qc_mode: str
    sample_count: int
    chromosome_count: int
    zero_row_jobs: int
    zero_row_autosomes: int
    warning_count: int
    error_count: int
    warnings: list[dict[str, str]]
    errors: list[dict[str, str]]
    sample_totals: list[dict[str, Any]]
    zero_row_autosome_records: list[dict[str, Any]]

    @property
    def is_valid(self) -> bool:
        return self.error_count == 0

    def to_dict(self) -> dict[str, Any]:
        data = asdict(self)
        data["is_valid"] = self.is_valid
        return data


def natural_chrom_key(chrom: Any) -> tuple[int, Any]:
    c = str(chrom).replace("chr", "")
    if c.isdigit():
        return (0, int(c))
    if c.upper() == "X":
        return (1, 23)
    if c.upper() == "Y":
        return (1, 24)
    if c.upper() in {"M", "MT"}:
        return (1, 25)
    return (2, c)


def issue(code: str, message: str, severity: str, sample_id: str = "", chromosome: str = "") -> dict[str, str]:
    out = {"code": code, "severity": severity, "message": message}
    if sample_id:
        out["sample_id"] = sample_id
    if chromosome:
        out["chromosome"] = chromosome
    return out


def run_build_qc(
    manifest: Path | str,
    out_dir: Path | str,
    *,
    qc_mode: str = "warn",
    print_summary: bool = True,
) -> BuildQcResult:
    manifest = Path(manifest)
    out_dir = Path(out_dir)
    if qc_mode not in QC_MODES:
        raise ValueError(f"Invalid qc_mode {qc_mode!r}; expected one of {sorted(QC_MODES)}")

    out_dir.mkdir(parents=True, exist_ok=True)
    with manifest.open(encoding="utf-8") as f:
        manifest_data = json.load(f)

    jobs = pd.DataFrame(manifest_data.get("jobs", []))
    merges = pd.DataFrame(manifest_data.get("sample_merges", []))

    if jobs.empty:
        warning = issue("NO_JOBS", "No chromosome jobs found in build manifest", "warning")
        result = BuildQcResult(
            manifest=str(manifest),
            out_dir=str(out_dir),
            qc_mode=qc_mode,
            sample_count=0,
            chromosome_count=0,
            zero_row_jobs=0,
            zero_row_autosomes=0,
            warning_count=1,
            error_count=0,
            warnings=[warning],
            errors=[],
            sample_totals=[],
            zero_row_autosome_records=[],
        )
        (out_dir / "build_qc_summary.json").write_text(json.dumps(result.to_dict(), indent=2), encoding="utf-8")
        return result

    jobs["flat_rows"] = pd.to_numeric(jobs["flat_rows"], errors="coerce").fillna(0).astype(int)
    jobs["chrom_order"] = jobs["chromosome"].map(natural_chrom_key)
    jobs = jobs.sort_values(["sample_id", "chrom_order"]).drop(columns=["chrom_order"])

    chrom_order = sorted(jobs["chromosome"].unique(), key=natural_chrom_key)
    row_matrix = jobs.pivot(index="sample_id", columns="chromosome", values="flat_rows")
    row_matrix = row_matrix[[c for c in chrom_order if c in row_matrix.columns]]

    skipped = jobs[jobs["flat_rows"].eq(0)].copy()
    autosome_zero = skipped[skipped["chromosome"].astype(str).str.fullmatch(r"[0-9]+")].copy()

    sample_totals = jobs.groupby("sample_id", as_index=False)["flat_rows"].sum().rename(
        columns={"flat_rows": "total_flat_rows"}
    )
    if not merges.empty and "variant_count" in merges.columns:
        cols = [c for c in ["sample_id", "variant_count", "merged_gds", "status"] if c in merges.columns]
        sample_totals = sample_totals.merge(merges[cols], on="sample_id", how="left")
        if "variant_count" in sample_totals.columns:
            sample_totals["flat_vs_merged_delta"] = (
                sample_totals["total_flat_rows"] - sample_totals["variant_count"].fillna(0).astype(int)
            )

    long_cols = [c for c in ["sample_id", "chromosome", "flat_rows", "status", "message", "summary_json"] if c in jobs.columns]
    long_counts = jobs[long_cols].copy()

    row_matrix.to_csv(out_dir / "flat_rows_by_sample_chromosome.csv")
    skipped.to_csv(out_dir / "skipped_chromosomes.csv", index=False)
    autosome_zero.to_csv(out_dir / "zero_row_autosomes.csv", index=False)
    sample_totals.to_csv(out_dir / "sample_total_rows.csv", index=False)
    long_counts.to_csv(out_dir / "all_job_counts_long.csv", index=False)

    warnings: list[dict[str, str]] = []
    errors: list[dict[str, str]] = []

    if qc_mode != "off":
        for _, row in autosome_zero.iterrows():
            sev = "error" if qc_mode == "strict" else "warning"
            rec = issue(
                "ZERO_ROW_AUTOSOME",
                "Sample/chromosome produced zero flat rows after carried-variant and GeneList filtering",
                sev,
                sample_id=str(row.get("sample_id", "")),
                chromosome=str(row.get("chromosome", "")),
            )
            if sev == "error":
                errors.append(rec)
            else:
                warnings.append(rec)

    zero_records_cols = [c for c in ["sample_id", "chromosome", "flat_rows", "status", "message", "summary_json"] if c in autosome_zero.columns]
    result = BuildQcResult(
        manifest=str(manifest),
        out_dir=str(out_dir),
        qc_mode=qc_mode,
        sample_count=int(jobs["sample_id"].nunique()),
        chromosome_count=int(len(chrom_order)),
        zero_row_jobs=int(len(skipped)),
        zero_row_autosomes=int(len(autosome_zero)),
        warning_count=len(warnings),
        error_count=len(errors),
        warnings=warnings,
        errors=errors,
        sample_totals=sample_totals.to_dict(orient="records"),
        zero_row_autosome_records=autosome_zero[zero_records_cols].to_dict(orient="records"),
    )

    (out_dir / "build_qc_summary.json").write_text(json.dumps(result.to_dict(), indent=2), encoding="utf-8")

    if print_summary:
        print("\nFlat row count matrix:")
        print(row_matrix.fillna("NA").to_string())
        print("\nSample total rows:")
        print(sample_totals.to_string(index=False))
        print("\nZero-row autosomes:")
        if autosome_zero.empty:
            print("None")
        else:
            print(autosome_zero[zero_records_cols].to_string(index=False))
        print(f"\nQC outputs written to: {out_dir}")
        if warnings:
            print(f"\nQC warnings: {len(warnings)}")
        if errors:
            print(f"\nQC errors: {len(errors)}")

    return result


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="QC build_manifest.json from IMPACT-SNV build-gds.")
    parser.add_argument("--manifest", required=True, type=Path, help="Path to build_manifest.json")
    parser.add_argument("--out-dir", required=True, type=Path, help="Directory for QC outputs")
    parser.add_argument("--qc-mode", choices=sorted(QC_MODES), default="warn")
    return parser


def run_qc_build(args: Any) -> int:
    result = run_build_qc(args.manifest, args.out_dir, qc_mode=getattr(args, "qc_mode", "warn"), print_summary=True)
    return 0 if result.is_valid else 2


def main(argv: Optional[list[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    return run_qc_build(args)


if __name__ == "__main__":
    raise SystemExit(main())
