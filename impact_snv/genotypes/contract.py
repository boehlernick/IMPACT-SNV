"""Genotype parquet layout contract for IMPACT-SNV.

This contract validates the extracted genotypes directory consumed by
`impact_snv.gds.flatten` and `impact_snv.gds.build`.
"""

from __future__ import annotations

from pathlib import Path
from typing import Iterable, List, Sequence

import pyarrow.parquet as pq


REQUIRED_BASE_COLUMNS: List[str] = ["chromosome", "position", "ref", "alt"]
DOSAGE_COLUMN_ALIASES: List[str] = ["dosages", "dosage"]


def chromosome_parquet_path(genotypes_dir: Path, chromosome: str) -> Path:
    return genotypes_dir / f"chromosome={chromosome}" / "data.parquet"


def validate_parquet_columns(columns: Iterable[str]) -> List[str]:
    names = set(columns)
    missing = [c for c in REQUIRED_BASE_COLUMNS if c not in names]
    missing.extend(c for c in DOSAGE_COLUMN_ALIASES if c not in names)
    return missing


def validate_dosage_vector_lengths(parquet_path: Path, sample_count: int) -> List[str]:
    if sample_count <= 0:
        return []

    table = pq.ParquetFile(parquet_path).read(columns=DOSAGE_COLUMN_ALIASES)
    errors: List[str] = []
    for column_name in DOSAGE_COLUMN_ALIASES:
        for row_index, value in enumerate(table.column(column_name).to_pylist(), start=1):
            if value is None:
                errors.append(
                    f"Chromosome dosage column {column_name} row {row_index} is null ({parquet_path})"
                )
                continue
            if len(value) != sample_count:
                errors.append(
                    f"Chromosome dosage column {column_name} row {row_index} has length {len(value)}; expected {sample_count} from samples.txt ({parquet_path})"
                )
                break
    return errors


def validate_genotypes_layout(genotypes_dir: Path, chromosomes: Sequence[str]) -> List[str]:
    """Return a list of contract violations for a genotype output directory."""
    errors: List[str] = []
    sample_count = 0
    samples_path = genotypes_dir / "samples.txt"
    if not samples_path.exists():
        errors.append(f"Missing required file: {samples_path}")
    else:
        sample_ids = [x.strip() for x in samples_path.read_text(encoding="utf-8").splitlines() if x.strip()]
        if not sample_ids:
            errors.append(f"samples.txt contains no sample IDs: {samples_path}")
        else:
            sample_count = len(sample_ids)

    for chrom in chromosomes:
        parquet_path = chromosome_parquet_path(genotypes_dir, str(chrom))
        if not parquet_path.exists():
            errors.append(f"Missing chromosome parquet: {parquet_path}")
            continue
        try:
            pf = pq.ParquetFile(parquet_path)
            missing = validate_parquet_columns(pf.schema_arrow.names)
            if missing:
                errors.append(
                    f"Chromosome {chrom} parquet is missing required columns: {missing} ({parquet_path})"
                )
                continue
            errors.extend(validate_dosage_vector_lengths(parquet_path, sample_count))
        except Exception as exc:
            errors.append(f"Unable to read chromosome parquet {parquet_path}: {exc}")

    return errors
