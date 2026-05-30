from pathlib import Path

import pandas as pd

from impact_snv.genotypes.contract import validate_genotypes_layout


def _write_genotype_parquet(
    path: Path,
    *,
    include_dosage: bool = True,
    vector_length: int = 1,
    rows: int = 1,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    records = []
    for idx in range(rows):
        row = {
            "chromosome": "1",
            "position": 101 + idx,
            "ref": "A",
            "alt": "G",
        }
        if include_dosage:
            dosage = [1.0] * vector_length
            row["dosages"] = dosage
            row["dosage"] = dosage
        records.append(row)
    columns = ["chromosome", "position", "ref", "alt"]
    if include_dosage:
        columns.extend(["dosages", "dosage"])
    pd.DataFrame(records, columns=columns).to_parquet(path, index=False, engine="pyarrow")


def test_validate_genotypes_layout_ok(tmp_path: Path) -> None:
    out = tmp_path / "genotypes"
    out.mkdir(parents=True, exist_ok=True)
    (out / "samples.txt").write_text("S001\n", encoding="utf-8")
    _write_genotype_parquet(out / "chromosome=1" / "data.parquet", include_dosage=True)

    errors = validate_genotypes_layout(out, ["1"])
    assert errors == []


def test_validate_genotypes_layout_reports_missing_dosage(tmp_path: Path) -> None:
    out = tmp_path / "genotypes"
    out.mkdir(parents=True, exist_ok=True)
    (out / "samples.txt").write_text("S001\n", encoding="utf-8")
    _write_genotype_parquet(out / "chromosome=1" / "data.parquet", include_dosage=False)

    errors = validate_genotypes_layout(out, ["1"])
    assert errors
    assert "dosages" in errors[0]
    assert "dosage" in errors[0]


def test_validate_genotypes_layout_allows_zero_row_chromosome(tmp_path: Path) -> None:
    out = tmp_path / "genotypes"
    out.mkdir(parents=True, exist_ok=True)
    (out / "samples.txt").write_text("S001\n", encoding="utf-8")
    _write_genotype_parquet(out / "chromosome=1" / "data.parquet", include_dosage=True, rows=0)

    errors = validate_genotypes_layout(out, ["1"])
    assert errors == []


def test_validate_genotypes_layout_reports_mismatched_dosage_vector_length(tmp_path: Path) -> None:
    out = tmp_path / "genotypes"
    out.mkdir(parents=True, exist_ok=True)
    (out / "samples.txt").write_text("S001\nS002\n", encoding="utf-8")
    _write_genotype_parquet(out / "chromosome=1" / "data.parquet", include_dosage=True, vector_length=1)

    errors = validate_genotypes_layout(out, ["1"])
    assert errors
    assert "expected 2 from samples.txt" in errors[0]
