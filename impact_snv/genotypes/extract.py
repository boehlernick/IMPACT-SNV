#!/usr/bin/env python3
"""Extract multi-sample VCF genotypes into IMPACT-SNV genotype parquet.

This module backs the production CLI command:

    impact-snv extract-genotypes

It writes the simple genotype parquet contract consumed by
impact_snv.gds.flatten/build:

    <out_dir>/samples.txt
    <out_dir>/chromosome=<chrom>/data.parquet
    <out_dir>/genotype_manifest.json

Each chromosome parquet contains one row per variant and list-valued dosage
vectors in the same sample order as samples.txt.

Canonical genotype parquet schema for v1.0.0:

    chromosome: string
    position: int64
    ref: string
    alt: string
    dosages: list<float32>
    dosage: list<float32>   # compatibility alias

The canonical column is `dosages`, matching the existing build-gds/flatten.py
contract. The singular `dosage` alias is retained temporarily for compatibility
with early adapter outputs and ad hoc inspection scripts.
"""

from __future__ import annotations

import json
import os
import re
import shutil
import subprocess
import sys
from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Optional, Sequence

import pyarrow as pa
import pyarrow.parquet as pq

SUPPORTED_REFERENCE_BUILDS = {"GRCh38"}
DEFAULT_CHROMS = [str(i) for i in range(1, 23)] + ["X", "Y"]
VERSION = "1.0.0a0"


@dataclass
class ChromosomeExtractResult:
    chromosome: str
    output_parquet: str
    row_count: int
    status: str
    message: str = ""


@dataclass
class GenotypeExtractResult:
    command: str
    version: str
    created_utc: str
    input_vcf: str
    input_index: Optional[str]
    out_dir: str
    reference_build: str
    sample_count: int
    sample_ids: list[str]
    chromosomes: list[str]
    chromosome_outputs: list[ChromosomeExtractResult]
    status: str
    warnings: list[dict[str, str]] = field(default_factory=list)
    errors: list[dict[str, str]] = field(default_factory=list)

    def to_dict(self) -> dict[str, Any]:
        data = asdict(self)
        data["chromosome_outputs"] = [asdict(x) for x in self.chromosome_outputs]
        return data


def utc_now_iso() -> str:
    return datetime.now(timezone.utc).isoformat()


def issue(code: str, message: str, severity: str = "warning") -> dict[str, str]:
    return {"code": code, "severity": severity, "message": message}


def resolve_executable(name: str, explicit: Optional[Path] = None) -> str:
    if explicit:
        p = Path(explicit).expanduser()
        if not p.exists() or not os.access(p, os.X_OK):
            raise FileNotFoundError(f"Executable for {name!r} is not executable: {p}")
        return str(p)
    found = shutil.which(name)
    if not found:
        raise FileNotFoundError(f"Required executable {name!r} was not found on PATH")
    return found


def existing_vcf_index(vcf: Path) -> Optional[Path]:
    for suffix in (".csi", ".tbi"):
        idx = Path(str(vcf) + suffix)
        if idx.exists():
            return idx.resolve()
    return None


def normalize_chrom(value: Any) -> str:
    text = str(value)
    if text.lower().startswith("chr"):
        text = text[3:]
    if text.endswith(".0"):
        text = text[:-2]
    return text


def parse_chromosomes(values: Optional[Sequence[str]]) -> list[str]:
    if not values:
        return list(DEFAULT_CHROMS)
    out: list[str] = []
    for value in values:
        value = str(value).strip()
        m = re.fullmatch(r"(\d+)-(\d+)", value)
        if m:
            start, end = int(m.group(1)), int(m.group(2))
            step = 1 if start <= end else -1
            out.extend(str(i) for i in range(start, end + step, step))
        else:
            out.append(normalize_chrom(value))
    seen: set[str] = set()
    final: list[str] = []
    for chrom in out:
        if chrom not in seen:
            seen.add(chrom)
            final.append(chrom)
    return final


def read_samples(vcf: Path, bcftools: str) -> list[str]:
    proc = subprocess.run([bcftools, "query", "-l", str(vcf)], text=True, capture_output=True)
    if proc.returncode != 0:
        raise RuntimeError(f"bcftools query -l failed for {vcf}:\n{proc.stderr}")
    samples = [x.strip() for x in proc.stdout.splitlines() if x.strip()]
    if not samples:
        raise RuntimeError(f"No sample IDs found in VCF: {vcf}")
    if len(samples) != len(set(samples)):
        raise RuntimeError("Duplicate sample IDs found in merged VCF; reheader before extract-genotypes")
    return samples


def allele_dosage_from_gt(gt: str) -> Optional[float]:
    """Convert GT to alternate-allele dosage.

    Missing alleles return None. Non-reference alleles are counted as 1 each, so:

        0/0 -> 0
        0/1, 1/0 -> 1
        1/1 -> 2
        haploid 0 -> 0
        haploid 1 -> 1
        1/2 -> 2

    This is deliberately conservative for multiallelic records; IMPACT-SNV
    merge/normalization should generally keep biallelic rows.
    """
    gt = (gt or "").strip()
    if not gt or gt in {".", "./.", ".|."}:
        return None
    gt = gt.split(":", 1)[0]
    alleles = re.split(r"[\/|]", gt)
    dosage = 0
    observed = False
    for allele in alleles:
        if allele in {"", "."}:
            return None
        try:
            val = int(allele)
        except ValueError:
            return None
        observed = True
        if val > 0:
            dosage += 1
    return float(dosage) if observed else None


def genotype_schema() -> pa.Schema:
    return pa.schema([
        ("chromosome", pa.string()),
        ("position", pa.int64()),
        ("ref", pa.string()),
        ("alt", pa.string()),
        ("dosages", pa.list_(pa.float32())),
        ("dosage", pa.list_(pa.float32())),
    ])


def write_empty_chromosome_parquet(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    schema = genotype_schema()
    pq.write_table(
        pa.Table.from_arrays([pa.array([], type=field.type) for field in schema], schema=schema),
        path,
    )


def extract_chromosome(
    *,
    vcf: Path,
    chrom: str,
    out_parquet: Path,
    sample_count: int,
    bcftools: str,
    progress_interval_records: int = 500_000,
) -> ChromosomeExtractResult:
    out_parquet.parent.mkdir(parents=True, exist_ok=True)
    fmt = "%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n"
    cmd = [bcftools, "query", "-r", chrom, "-f", fmt, str(vcf)]

    chromosomes: list[str] = []
    positions: list[int] = []
    refs: list[str] = []
    alts: list[str] = []
    dosages: list[list[Optional[float]]] = []
    row_count = 0

    proc = subprocess.Popen(cmd, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    assert proc.stdout is not None
    try:
        for line in proc.stdout:
            line = line.rstrip("\n")
            if not line:
                continue
            fields = line.split("\t")
            if len(fields) < 4 + sample_count:
                raise RuntimeError(
                    f"Malformed bcftools query row for chromosome {chrom}: "
                    f"expected at least {4 + sample_count} fields, got {len(fields)}"
                )
            chromosomes.append(normalize_chrom(fields[0]))
            positions.append(int(fields[1]))
            refs.append(fields[2])
            alts.append(fields[3])
            dosages.append([allele_dosage_from_gt(x) for x in fields[4:4 + sample_count]])
            row_count += 1
            if progress_interval_records and row_count % progress_interval_records == 0:
                print(f"  {chrom}: extracted {row_count:,} genotype rows", flush=True)
    finally:
        proc.stdout.close()

    stderr = proc.stderr.read() if proc.stderr else ""
    rc = proc.wait()
    if proc.stderr:
        proc.stderr.close()
    if rc != 0:
        return ChromosomeExtractResult(chrom, str(out_parquet), 0, "error", stderr.strip())

    if row_count == 0:
        write_empty_chromosome_parquet(out_parquet)
        return ChromosomeExtractResult(chrom, str(out_parquet), 0, "ok", "No variants found for chromosome")

    dosage_array = pa.array(dosages, type=pa.list_(pa.float32()))
    table = pa.table({
        "chromosome": pa.array(chromosomes, type=pa.string()),
        "position": pa.array(positions, type=pa.int64()),
        "ref": pa.array(refs, type=pa.string()),
        "alt": pa.array(alts, type=pa.string()),
        "dosages": dosage_array,
        "dosage": dosage_array,
    }, schema=genotype_schema())
    pq.write_table(table, out_parquet)
    return ChromosomeExtractResult(chrom, str(out_parquet), row_count, "ok")


def write_manifest(result: GenotypeExtractResult, manifest: Path) -> None:
    manifest.parent.mkdir(parents=True, exist_ok=True)
    manifest.write_text(json.dumps(result.to_dict(), indent=2, sort_keys=True), encoding="utf-8")


def run_extract_genotypes(args: Any) -> int:
    reference_build = getattr(args, "reference_build", "GRCh38")
    if reference_build not in SUPPORTED_REFERENCE_BUILDS:
        print(
            f"ERROR: IMPACT-SNV v1.0.0 supports only --reference-build GRCh38; got {reference_build}",
            file=sys.stderr,
        )
        return 2

    input_vcf = Path(args.input_vcf).expanduser().resolve()
    if not input_vcf.exists() or not input_vcf.is_file():
        print(f"ERROR: input VCF not found: {input_vcf}", file=sys.stderr)
        return 2
    input_index = existing_vcf_index(input_vcf)
    if input_index is None:
        print(f"ERROR: input VCF index not found; expected {input_vcf}.csi or {input_vcf}.tbi", file=sys.stderr)
        return 2

    out_dir = Path(args.out_dir).expanduser().resolve()
    force = bool(getattr(args, "force", False))
    if out_dir.exists() and any(out_dir.iterdir()) and not force:
        print(f"ERROR: output directory exists and is not empty; use --force: {out_dir}", file=sys.stderr)
        return 2
    if out_dir.exists() and force:
        shutil.rmtree(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    bcftools = resolve_executable("bcftools", Path(args.bcftools) if getattr(args, "bcftools", None) else None)
    chromosomes = parse_chromosomes(getattr(args, "chromosomes", None))
    samples = read_samples(input_vcf, bcftools)
    (out_dir / "samples.txt").write_text("\n".join(samples) + "\n", encoding="utf-8")

    print(f"Extracting genotypes from {input_vcf}", flush=True)
    print(f"  samples: {len(samples)}", flush=True)
    print(f"  chromosomes: {' '.join(chromosomes)}", flush=True)
    print(f"  out_dir: {out_dir}", flush=True)

    chrom_results: list[ChromosomeExtractResult] = []
    for chrom in chromosomes:
        out_parquet = out_dir / f"chromosome={chrom}" / "data.parquet"
        print(f"Extracting chromosome {chrom}...", flush=True)
        res = extract_chromosome(
            vcf=input_vcf,
            chrom=chrom,
            out_parquet=out_parquet,
            sample_count=len(samples),
            bcftools=bcftools,
            progress_interval_records=int(getattr(args, "progress_interval_records", 500_000) or 500_000),
        )
        chrom_results.append(res)
        print(f"  {chrom}: {res.status}, rows={res.row_count:,}", flush=True)
        if res.status == "error" and getattr(args, "qc_mode", "warn") == "strict":
            break

    errors = [
        issue("CHROMOSOME_EXTRACT_ERROR", f"{r.chromosome}: {r.message}", "error")
        for r in chrom_results
        if r.status == "error"
    ]
    warnings = [
        issue("ZERO_ROW_CHROMOSOME", f"Chromosome {r.chromosome} produced zero genotype rows")
        for r in chrom_results
        if r.status == "ok" and r.row_count == 0
    ]
    status = "ok" if not errors else "failed"

    result = GenotypeExtractResult(
        command="impact-snv extract-genotypes",
        version=VERSION,
        created_utc=utc_now_iso(),
        input_vcf=str(input_vcf),
        input_index=str(input_index) if input_index else None,
        out_dir=str(out_dir),
        reference_build=reference_build,
        sample_count=len(samples),
        sample_ids=samples,
        chromosomes=chromosomes,
        chromosome_outputs=chrom_results,
        status=status,
        warnings=warnings,
        errors=errors,
    )
    manifest = Path(getattr(args, "manifest_json", None) or (out_dir / "genotype_manifest.json"))
    write_manifest(result, manifest)
    print(f"Genotype manifest written: {manifest}", flush=True)

    if status != "ok":
        print(f"Genotype extraction failed with {len(errors)} error(s)", file=sys.stderr)
        return 1
    print("Genotype extraction completed successfully.", flush=True)
    return 0


__all__ = ["run_extract_genotypes", "GenotypeExtractResult", "ChromosomeExtractResult"]
