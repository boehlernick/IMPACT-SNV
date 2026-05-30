from __future__ import annotations

import json
import shutil
import subprocess
import sys
from functools import lru_cache
from pathlib import Path

import pytest

from impact_snv.cli import main


REPO_ROOT = Path(__file__).resolve().parents[1]
GENEBREAKER_ROOT = REPO_ROOT / "tests" / "genebreaker_vcfs" / "merged_vcf_test"


@lru_cache(maxsize=1)
def _r_gds_runtime() -> tuple[bool, str]:
    rscript = shutil.which("Rscript")
    if not rscript:
        return False, "Rscript not available"

    proc = subprocess.run(
        [
            rscript,
            "-e",
            "suppressPackageStartupMessages({library(SeqArray); library(gdsfmt)}); cat('R_GDS_OK\\n')",
        ],
        text=True,
        capture_output=True,
    )
    if proc.returncode != 0:
        detail = (proc.stderr or proc.stdout or "").strip()
        return False, detail or "SeqArray/gdsfmt runtime unavailable"
    return True, rscript


def _require_r_gds_runtime() -> str:
    ok, detail = _r_gds_runtime()
    if not ok:
        pytest.skip(f"requires Rscript with SeqArray and gdsfmt: {detail}")
    return detail


def _fixture_path(*parts: str) -> Path:
    path = GENEBREAKER_ROOT.joinpath(*parts)
    if not path.exists():
        pytest.skip(f"required release fixture is missing: {path}")
    return path


def test_release_cli_end_to_end_build_finalize_validate_and_qc(tmp_path: Path) -> None:
    rscript = _require_r_gds_runtime()

    annotated_dir = _fixture_path("merged_case1_case2_output.annotated")
    genotypes_dir = _fixture_path("merged_case1_case2_output.genotypes")
    gene_list = _fixture_path("GeneList.txt")

    sample_id = "Case1_mother"
    build_out = tmp_path / "build"
    finalize_out = tmp_path / "final"
    qc_recheck_out = tmp_path / "qc_recheck"

    rc = main(
        [
            "build-gds",
            "--annotated-dir",
            str(annotated_dir),
            "--genotypes-dir",
            str(genotypes_dir),
            "--gene-list",
            str(gene_list),
            "--out-dir",
            str(build_out),
            "--samples",
            sample_id,
            "--chromosomes",
            "1",
            "X",
            "--python",
            sys.executable,
            "--rscript",
            rscript,
            "--qc-mode",
            "strict",
            "--no-optimize",
        ]
    )
    assert rc == 0

    build_manifest_path = build_out / "build_manifest.json"
    assert build_manifest_path.exists()
    build_manifest = json.loads(build_manifest_path.read_text(encoding="utf-8"))
    assert build_manifest["command"] == "impact-snv build-gds"
    assert build_manifest["samples"] == [sample_id]
    assert build_manifest["chromosomes"] == ["1", "X"]
    assert len(build_manifest["sample_merges"]) == 1
    assert build_manifest["sample_merges"][0]["status"] == "ok"
    assert build_manifest["sample_merges"][0]["variant_count"] > 0

    preprior_gds = build_out / "gds_merged" / f"{sample_id}_SNV_IMPACT.preprioritization.gds"
    assert preprior_gds.exists()

    build_qc_summary_path = build_out / "qc" / "build_qc_summary.json"
    assert build_qc_summary_path.exists()
    build_qc_summary = json.loads(build_qc_summary_path.read_text(encoding="utf-8"))
    assert build_qc_summary["is_valid"] is True
    assert build_qc_summary["sample_count"] == 1
    assert build_qc_summary["warning_count"] == 0
    assert build_qc_summary["error_count"] == 0
    assert build_qc_summary["zero_row_autosomes"] == 0

    rc = main(
        [
            "finalize-gds",
            "--input-dir",
            str(build_out / "gds_merged"),
            "--gene-list",
            str(gene_list),
            "--out-dir",
            str(finalize_out),
            "--samples",
            sample_id,
            "--rscript",
            rscript,
            "--qc-mode",
            "strict",
            "--no-optimize",
        ]
    )
    assert rc == 0

    final_gds = finalize_out / f"{sample_id}_SNV_IMPACT.gds"
    assert final_gds.exists()

    finalize_manifest_path = finalize_out / "finalize_manifest.json"
    assert finalize_manifest_path.exists()
    finalize_manifest = json.loads(finalize_manifest_path.read_text(encoding="utf-8"))
    assert finalize_manifest["ok_count"] == 1
    assert finalize_manifest["failed_count"] == 0
    result = finalize_manifest["results"][0]
    assert result["status"] == "ok"
    assert result["validation"]["is_valid"] is True
    assert result["validation"]["missing_nodes"] == []
    assert result["validation"]["variant_count"] > 0
    assert result["validation"]["varinfo_preview"]

    rc = main(
        [
            "validate-gds",
            "--gds",
            str(final_gds),
            "--rscript",
            rscript,
            "--qc-mode",
            "strict",
        ]
    )
    assert rc == 0

    rc = main(
        [
            "qc-build",
            "--manifest",
            str(build_manifest_path),
            "--out-dir",
            str(qc_recheck_out),
            "--qc-mode",
            "strict",
        ]
    )
    assert rc == 0

    qc_recheck_summary_path = qc_recheck_out / "build_qc_summary.json"
    assert qc_recheck_summary_path.exists()
    qc_recheck_summary = json.loads(qc_recheck_summary_path.read_text(encoding="utf-8"))
    assert qc_recheck_summary["is_valid"] is True
    assert qc_recheck_summary["sample_count"] == 1
    assert qc_recheck_summary["warning_count"] == 0
    assert qc_recheck_summary["error_count"] == 0