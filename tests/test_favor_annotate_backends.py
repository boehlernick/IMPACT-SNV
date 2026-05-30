from argparse import Namespace
from pathlib import Path
import json

import pandas as pd

from impact_snv.favor.annotate import run_favor_annotate


def _write_parquet(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame([
        {
            "chromosome": "1",
            "position": 101,
            "ref": "A",
            "alt": "G",
            "clnsig": "Pathogenic",
        }
    ]).to_parquet(path, index=False, engine="pyarrow")


def test_legacy_backend_stages_and_writes_manifest(tmp_path: Path) -> None:
    legacy_annotated = tmp_path / "legacy.annotated"
    legacy_genotypes = tmp_path / "legacy.genotypes"
    out_dir = tmp_path / "out"

    _write_parquet(legacy_annotated / "chromosome=1" / "data.parquet")
    _write_parquet(legacy_genotypes / "chromosome=1" / "data.parquet")
    (legacy_genotypes / "samples.txt").write_text("S001\n", encoding="utf-8")

    args = Namespace(
        backend="legacy-favorannotator",
        out_dir=out_dir,
        out_prefix="caseA",
        reference_build="GRCh38",
        legacy_annotated_dir=legacy_annotated,
        legacy_genotypes_dir=legacy_genotypes,
        legacy_stage_mode="symlink",
        manifest_json=None,
        force=False,
    )

    rc = run_favor_annotate(args)
    assert rc == 0

    staged_annotated = out_dir / "caseA.annotated"
    staged_genotypes = out_dir / "caseA.genotypes"
    manifest = out_dir / "caseA.annotate_manifest.json"

    assert staged_annotated.exists()
    assert staged_genotypes.exists()
    assert manifest.exists()


def test_skeleton_backend_writes_manifest_with_dry_run(tmp_path: Path, monkeypatch) -> None:
    input_gds = tmp_path / "input.gds"
    input_gds.write_text("placeholder", encoding="utf-8")
    out_dir = tmp_path / "out"

    def _fake_run(*args, **kwargs):
        env = kwargs.get("env", {})
        out = Path(env["FAVORCLI_OUTPUT_DIR"])
        out.mkdir(parents=True, exist_ok=True)
        (out / "variant_identity.tsv").write_text("variant.id\n1\n", encoding="utf-8")
        (out / "favorcli_validation_report.txt").write_text("ok\n", encoding="utf-8")
        (out / "favorcli_run_manifest.json").write_text("{}\n", encoding="utf-8")
        return None

    monkeypatch.setattr("impact_snv.favor.annotate.subprocess.run", _fake_run)

    args = Namespace(
        backend="favor-cli-skeleton",
        out_dir=out_dir,
        out_prefix="caseB",
        reference_build="GRCh38",
        skeleton_input_file=input_gds,
        skeleton_input_type="gds",
        skeleton_dry_run=True,
        favor_bin="favor",
        favor_database_file=None,
        favor_database_version=None,
        skeleton_threads=4,
        skeleton_memory_budget_gb=8,
        manifest_json=None,
        force=False,
    )

    rc = run_favor_annotate(args)
    assert rc == 0

    manifest = out_dir / "caseB.annotate_manifest.json"
    assert manifest.exists()
    payload = json.loads(manifest.read_text(encoding="utf-8"))
    assert payload["backend"] == "favor-cli-skeleton"


def test_favor_cli_backend_allows_annotation_only_success_and_records_warnings(tmp_path: Path, monkeypatch) -> None:
    ingested_dir = tmp_path / "caseC.ingested"
    ingested_dir.mkdir(parents=True, exist_ok=True)
    out_dir = tmp_path / "out"
    favor_bin = tmp_path / "favor"
    favor_bin.write_text("#!/bin/sh\nexit 0\n", encoding="utf-8")

    def _fake_run_logged_command_with_progress(command, **kwargs):
        annotated_dir = Path(kwargs["watch_paths"][0])
        _write_parquet(annotated_dir / "chromosome=1" / "data.parquet")
        return 0

    monkeypatch.setattr(
        "impact_snv.favor.annotate.run_logged_command_with_progress",
        _fake_run_logged_command_with_progress,
    )

    args = Namespace(
        backend="favor-cli",
        ingested_dir=ingested_dir,
        out_dir=out_dir,
        out_prefix="caseC",
        reference_build="GRCh38",
        favor_bin=str(favor_bin),
        force=False,
        manifest_json=None,
        progress_interval_seconds=60,
        quiet=False,
        tail_log_lines=0,
        max_progress_log_line_chars=300,
        progress_mode="normal",
    )

    rc = run_favor_annotate(args)
    assert rc == 0

    manifest = out_dir / "caseC.annotate_manifest.json"
    assert manifest.exists()
    payload = json.loads(manifest.read_text(encoding="utf-8"))
    warning_codes = {warning["code"] for warning in payload["warnings"]}

    assert payload["backend"] == "favor-cli"
    assert payload["status"] == "ok"
    assert payload["annotated_dir"].endswith("caseC.annotated")
    assert payload["genotypes_dir"].endswith("caseC.genotypes")
    assert payload["samples_file"].endswith("caseC.genotypes/samples.txt")
    assert payload["annotated_chromosome_partitions"]
    assert payload["genotype_chromosome_partitions"] == []
    assert payload["favor_command"][1] == "annotate"
    assert payload["stdout_log"].endswith("caseC.favor_annotate.stdout.log")
    assert payload["stderr_log"].endswith("caseC.favor_annotate.stderr.log")
    assert warning_codes == {"GENOTYPES_OUTPUT_NOT_CREATED", "GENOTYPE_SAMPLES_FILE_NOT_CREATED"}
