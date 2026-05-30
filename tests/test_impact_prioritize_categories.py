from __future__ import annotations

import os
import shutil
import subprocess
from functools import lru_cache
from pathlib import Path

import pytest


REPO_ROOT = Path(__file__).resolve().parents[1]
PRIORITIZE_SCRIPT = REPO_ROOT / "impact_snv" / "resources" / "impact_prioritize_gds.R"


@lru_cache(maxsize=1)
def _rscript_path() -> str:
    rscript = shutil.which("Rscript")
    if not rscript:
        pytest.skip("Rscript not available")
    return rscript


def test_exonic_category_classifier_is_token_safe() -> None:
    rscript = _rscript_path()
    program = f"""
source('{PRIORITIZE_SCRIPT.as_posix()}')
cases <- c(
  'nonframeshift',
  'frameshift',
  'nonsynonymous',
  'synonymous',
  'stoploss',
  'stopgain',
  ' InFrame_Insertion | STOP_LOST '
)
for (value in cases) {{
  classes <- paste(sort(classify_exonic_category(value)), collapse = ';')
  tier2 <- category_has_class(list(value), TIER2_EXONIC_CLASSES)
  tier3 <- category_has_class(list(value), TIER3_EXONIC_CLASSES)
  cat(value, '\t', classes, '\t', tier2, '\t', tier3, '\n', sep = '')
}}
"""
    proc = subprocess.run(
        [rscript, "-"],
        input=program,
        text=True,
        capture_output=True,
        env={**os.environ, "RENV_CONFIG_AUTOLOADER_ENABLED": "FALSE"},
    )
    if proc.returncode != 0:
        pytest.fail(proc.stderr or proc.stdout or "R category classifier test failed")

    rows: dict[str, tuple[str, str, str]] = {}
    for line in proc.stdout.strip().splitlines():
        if "\t" not in line:
            continue
        value, classes, tier2, tier3 = line.split("\t")
        rows[value] = (classes, tier2, tier3)

    assert rows["nonframeshift"] == ("nonframeshift", "FALSE", "TRUE")
    assert rows["frameshift"] == ("frameshift", "TRUE", "FALSE")
    assert rows["nonsynonymous"] == ("nonsynonymous", "FALSE", "TRUE")
    assert rows["synonymous"] == ("synonymous", "FALSE", "FALSE")
    assert rows["stoploss"] == ("stoploss", "FALSE", "TRUE")
    assert rows["stopgain"] == ("stopgain", "TRUE", "FALSE")
    assert rows[" InFrame_Insertion | STOP_LOST "] == (
        "nonframeshift_insertion;stoploss",
        "FALSE",
        "TRUE",
    )