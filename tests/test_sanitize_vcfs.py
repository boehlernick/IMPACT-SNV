from __future__ import annotations

import csv
import json
import shutil
import subprocess
from pathlib import Path

import pytest

from impact_snv.cli import main
from impact_snv.merge.sanitize import assign_unique_prefixes, build_sanitized_sample_ids, derive_sample_prefix


def _bcftools_runtime() -> tuple[str, str]:
    bcftools = shutil.which("bcftools")
    bgzip = shutil.which("bgzip")
    if not bcftools or not bgzip:
        pytest.skip("requires bcftools and bgzip")
    return bcftools, bgzip


def _write_minimal_trio_vcf(path: Path, pos: int) -> None:
    path.write_text(
        "\n".join(
            [
                "##fileformat=VCFv4.2",
                "##contig=<ID=1,length=1000>",
                "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tproband\tmother\tfather",
                f"1\t{pos}\t.\tA\tG\t.\tPASS\t.\tGT\t0/1\t0/0\t0/0",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


def test_derive_sample_prefix_prefers_case_id() -> None:
    assert derive_sample_prefix(Path("Case1_WAS_Merged_GeneBreaker (1).vcf.gz")) == "Case1"
    assert derive_sample_prefix(Path("run_2026_panel.vcf.gz")) == "run_2026_panel"


def test_build_sanitized_sample_ids_handles_local_and_global_duplicates() -> None:
    seen: set[str] = set()
    first = build_sanitized_sample_ids(["proband", "mother", "father"], "Case1", seen)
    second = build_sanitized_sample_ids(["proband", "proband"], "Case1_2", seen)
    assert first == ("Case1_proband", "Case1_mother", "Case1_father")
    assert second == ("Case1_2_proband", "Case1_2_proband_2")


def test_assign_unique_prefixes_disambiguates_duplicate_case_tokens() -> None:
    vcfs = [
        Path("Case1_alpha.vcf.gz"),
        Path("Case1_beta.vcf.gz"),
        Path("Case2_gamma.vcf.gz"),
    ]
    prefixes = assign_unique_prefixes(vcfs)
    assert prefixes[vcfs[0]] == "Case1"
    assert prefixes[vcfs[1]] == "Case1_2"
    assert prefixes[vcfs[2]] == "Case2"


def test_sanitize_vcfs_reheaders_duplicate_trio_inputs(tmp_path: Path) -> None:
    bcftools, bgzip = _bcftools_runtime()

    case1 = tmp_path / "Case1_alpha (1).vcf"
    case2 = tmp_path / "Case2_beta.vcf"
    _write_minimal_trio_vcf(case1, 101)
    _write_minimal_trio_vcf(case2, 202)

    out_dir = tmp_path / "sanitized"
    rc = main(
        [
            "sanitize-vcfs",
            "--vcfs",
            str(case1),
            str(case2),
            "--out-dir",
            str(out_dir),
            "--bcftools",
            bcftools,
            "--bgzip",
            bgzip,
        ]
    )
    assert rc == 0

    manifest_tsv = out_dir / "sanitized_vcfs.tsv"
    manifest_json = out_dir / "sanitize_manifest.json"
    assert manifest_tsv.exists()
    assert manifest_json.exists()

    manifest = json.loads(manifest_json.read_text(encoding="utf-8"))
    assert manifest["command"] == "impact-snv sanitize-vcfs"
    assert manifest["sample_count"] == 6
    assert len(manifest["records"]) == 2

    rows = list(csv.DictReader(manifest_tsv.read_text(encoding="utf-8").splitlines(), delimiter="\t"))
    assert len(rows) == 2
    out_paths = [Path(row["vcf_path"]) for row in rows]
    assert all(path.exists() for path in out_paths)

    expected = {
        "Case1.sanitized.vcf.gz": ["Case1_proband", "Case1_mother", "Case1_father"],
        "Case2.sanitized.vcf.gz": ["Case2_proband", "Case2_mother", "Case2_father"],
    }
    for out_vcf in out_paths:
        proc = subprocess.run([bcftools, "query", "-l", str(out_vcf)], text=True, capture_output=True, check=True)
        assert proc.stdout.splitlines() == expected[out_vcf.name]