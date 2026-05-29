#!/usr/bin/env python3
"""Common FAVOR CLI wrapper utilities for IMPACT-SNV.

This module intentionally keeps terminal progress compact by default while
preserving complete stdout/stderr logs on disk. Long-running FAVOR commands are
run with subprocess.Popen so users can see periodic heartbeat updates.
"""

from __future__ import annotations

import json
import re
import shutil
import subprocess
import threading
import time
from collections import deque
from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Optional, Sequence

SUPPORTED_REFERENCE_BUILDS = {"GRCh38"}
NO_ACTIVITY_WARNING_SECONDS = 30 * 60
NOISY_TERMINAL_LOG_SUBSTRINGS = (
    "Plan:",
    "ProjectionExec:",
    "HashJoinExec:",
    "RepartitionExec:",
    "DataSourceExec:",
    "TableScan:",
)
CHROM_STAGE_RE = re.compile(r"\bchr([0-9XYM]+)\s+\((\d+)/24\)", re.IGNORECASE)


@dataclass
class FavorRunResult:
    command: str
    favor_command: list[str]
    cwd: str
    out_dir: str
    out_prefix: str
    reference_build: str
    favor_bin: str
    force: bool
    return_code: int
    start_time_utc: str
    end_time_utc: str
    stdout_log: str
    stderr_log: str
    progress_interval_seconds: int
    quiet: bool
    tail_log_lines: int
    max_progress_log_line_chars: int
    progress_mode: str
    removed_paths: list[str] = field(default_factory=list)
    observed_output_sizes_bytes: dict[str, int] = field(default_factory=dict)
    observed_output_latest_mtime_utc: dict[str, Optional[str]] = field(default_factory=dict)
    status: str = "ok"
    message: str = ""

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


def utc_now_iso() -> str:
    return datetime.now(timezone.utc).isoformat()


def seconds_to_hhmmss(seconds: float) -> str:
    seconds_i = max(0, int(seconds))
    return f"{seconds_i // 3600:02d}:{(seconds_i % 3600) // 60:02d}:{seconds_i % 60:02d}"


def seconds_to_short(seconds: Optional[float]) -> str:
    if seconds is None:
        return "none"
    seconds_i = max(0, int(seconds))
    if seconds_i < 60:
        return f"{seconds_i}s ago"
    if seconds_i < 3600:
        return f"{seconds_i // 60}m {seconds_i % 60}s ago"
    return f"{seconds_i // 3600}h {(seconds_i % 3600) // 60}m ago"


def mtime_to_utc(mtime: Optional[float]) -> Optional[str]:
    if mtime is None:
        return None
    return datetime.fromtimestamp(mtime, tz=timezone.utc).isoformat()


def validate_reference_build(reference_build: str) -> None:
    if reference_build not in SUPPORTED_REFERENCE_BUILDS:
        raise ValueError(
            f"IMPACT-SNV v1.0.0 supports only reference-build GRCh38; got {reference_build!r}"
        )


def resolve_favor_bin(favor_bin: str) -> str:
    if Path(favor_bin).parent != Path("."):
        p = Path(favor_bin).expanduser()
        if not p.exists() or not p.is_file():
            raise FileNotFoundError(f"FAVOR executable not found: {p}")
        return str(p)
    if not shutil.which(favor_bin):
        raise FileNotFoundError(f"FAVOR executable {favor_bin!r} was not found on PATH")
    return favor_bin


def existing_vcf_index(input_vcf: Path) -> Optional[Path]:
    for suffix in (".csi", ".tbi"):
        idx = Path(str(input_vcf) + suffix)
        if idx.exists():
            return idx.resolve()
    return None


def infer_prefix_from_vcf(input_vcf: Path) -> str:
    name = input_vcf.name
    for suffix in (".vcf.gz", ".vcf.bgz", ".vcf", ".bcf"):
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return input_vcf.stem


def ensure_out_dir(out_dir: Path) -> Path:
    out_dir = Path(out_dir).expanduser().resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "logs").mkdir(parents=True, exist_ok=True)
    return out_dir


def remove_existing_outputs(paths: Sequence[Path], *, force: bool) -> list[str]:
    existing = [p for p in paths if p.exists()]
    if existing and not force:
        raise FileExistsError(
            "Existing FAVOR output path(s); use --force to remove before running:\n"
            + "\n".join(str(p) for p in existing)
        )
    removed: list[str] = []
    for p in existing:
        if p.is_dir():
            shutil.rmtree(p)
        else:
            p.unlink()
        removed.append(str(p.resolve()))
    return removed


def directory_size_bytes(path: Path) -> int:
    if not path.exists():
        return 0
    if path.is_file():
        try:
            return path.stat().st_size
        except FileNotFoundError:
            return 0
    total = 0
    for p in path.rglob("*"):
        try:
            if p.is_file():
                total += p.stat().st_size
        except FileNotFoundError:
            pass
    return total


def latest_mtime(path: Path) -> Optional[float]:
    if not path.exists():
        return None
    mtimes: list[float] = []
    try:
        mtimes.append(path.stat().st_mtime)
    except FileNotFoundError:
        return None
    if path.is_dir():
        for p in path.rglob("*"):
            try:
                mtimes.append(p.stat().st_mtime)
            except FileNotFoundError:
                pass
    return max(mtimes) if mtimes else None


def human_size(num_bytes: int) -> str:
    value = float(num_bytes)
    for unit in ("B", "KB", "MB", "GB", "TB"):
        if value < 1024.0 or unit == "TB":
            if unit == "B":
                return f"{int(value)} {unit}"
            return f"{value:.1f} {unit}"
        value /= 1024.0
    return f"{num_bytes} B"


def chromosome_partition_dirs(base: Path) -> list[str]:
    if not base.exists() or not base.is_dir():
        return []
    return sorted(str(p.resolve()) for p in base.iterdir() if p.is_dir() and p.name.startswith("chromosome="))


def tail_file(path: Path, n: int, *, max_scan_lines: int = 500) -> list[str]:
    """Read the last N lines without loading very large logs entirely."""
    if n <= 0 or not path.exists() or not path.is_file():
        return []
    try:
        with path.open("r", encoding="utf-8", errors="replace") as fh:
            return list(deque(fh, maxlen=max(n, min(max_scan_lines, n))))[-n:]
    except Exception:
        return []


def recent_lines(path: Path, *, n: int = 200) -> list[str]:
    return tail_file(path, n)


def sanitize_terminal_log_line(line: str, *, max_chars: int) -> Optional[str]:
    stripped = line.strip()
    if not stripped:
        return None
    if any(token in stripped for token in NOISY_TERMINAL_LOG_SUBSTRINGS):
        return None
    max_chars = max(40, int(max_chars or 300))
    if len(stripped) > max_chars:
        return stripped[:max_chars] + "... [truncated; see log for full line]"
    return stripped


def useful_tail_lines(path: Path, *, n: int, max_chars: int) -> list[str]:
    if n <= 0:
        return []
    # Scan more lines than requested because noisy DataFusion plan lines are skipped.
    candidates = recent_lines(path, n=max(100, n * 20))
    useful: list[str] = []
    for line in candidates:
        cleaned = sanitize_terminal_log_line(line, max_chars=max_chars)
        if cleaned is not None:
            useful.append(cleaned)
    return useful[-n:]


def current_chromosome_stage(stderr_log: Path) -> Optional[str]:
    for line in reversed(recent_lines(stderr_log, n=300)):
        m = CHROM_STAGE_RE.search(line)
        if m:
            chrom = m.group(1).upper()
            if chrom == "M":
                chrom = "M"
            return f"chr{chrom} ({m.group(2)}/24)"
    return None


def summarize_watch_paths(watch_paths: Sequence[Path], *, now: Optional[float] = None) -> dict[str, Any]:
    now = time.time() if now is None else now
    summary: dict[str, Any] = {
        "annotated_parts": None,
        "genotypes_parts": None,
        "ingested_exists": None,
        "latest_mtime": None,
        "latest_activity_ago": None,
        "paths": [],
    }
    latest_seen: Optional[float] = None
    for path in watch_paths:
        exists = path.exists()
        parts = len(chromosome_partition_dirs(path)) if exists and path.is_dir() else 0
        size = directory_size_bytes(path)
        mt = latest_mtime(path)
        if mt is not None:
            latest_seen = mt if latest_seen is None else max(latest_seen, mt)
        lower = path.name.lower()
        if lower.endswith(".annotated"):
            summary["annotated_parts"] = parts if exists else "pending"
        elif lower.endswith(".genotypes"):
            summary["genotypes_parts"] = parts if exists else "pending"
        elif lower.endswith(".ingested"):
            summary["ingested_exists"] = exists
        summary["paths"].append({"path": path, "exists": exists, "parts": parts, "size": size, "mtime": mt})
    summary["latest_mtime"] = latest_seen
    summary["latest_activity_ago"] = None if latest_seen is None else now - latest_seen
    return summary


def print_compact_progress(
    *,
    progress_label: str,
    stderr_log: Path,
    start_time: float,
    watch_paths: Sequence[Path],
) -> None:
    now = time.time()
    elapsed = seconds_to_hhmmss(now - start_time)
    stage = current_chromosome_stage(stderr_log)
    summary = summarize_watch_paths(watch_paths, now=now)
    fields = [f"[{progress_label}] running {elapsed}"]
    if stage:
        fields.append(stage)
    if summary.get("annotated_parts") is not None:
        fields.append(f"annotated parts: {summary['annotated_parts']}")
    if summary.get("genotypes_parts") is not None:
        gt = summary["genotypes_parts"]
        fields.append("genotypes: pending" if gt == "pending" else f"genotype parts: {gt}")
    if summary.get("ingested_exists") is not None:
        fields.append(f"ingested: {'yes' if summary['ingested_exists'] else 'pending'}")
    fields.append(f"activity: {seconds_to_short(summary.get('latest_activity_ago'))}")
    print(" | ".join(fields), flush=True)


def print_expanded_progress(
    *,
    progress_label: str,
    stderr_log: Path,
    stdout_log: Path,
    start_time: float,
    watch_paths: Sequence[Path],
    tail_log_lines: int,
    max_progress_log_line_chars: int,
    show_paths: bool,
    command: Sequence[str],
    cwd: Path,
) -> None:
    now = time.time()
    print(f"[{progress_label}] still running after {seconds_to_hhmmss(now - start_time)}", flush=True)
    stage = current_chromosome_stage(stderr_log)
    if stage:
        print(f"  current stage: {stage}", flush=True)
    if show_paths:
        print(f"  command: {' '.join(map(str, command))}", flush=True)
        print(f"  cwd: {cwd}", flush=True)
    for item in summarize_watch_paths(watch_paths, now=now)["paths"]:
        path = item["path"]
        exists = item["exists"]
        if not exists:
            print(f"  {path.name}: not created yet", flush=True)
        else:
            part_word = "partition" if item["parts"] == 1 else "partitions"
            print(
                f"  {path.name}: exists, {item['parts']} chromosome {part_word}, {human_size(item['size'])}",
                flush=True,
            )
    latest_ago = summarize_watch_paths(watch_paths, now=now).get("latest_activity_ago")
    print(f"  latest filesystem activity: {seconds_to_short(latest_ago)}", flush=True)
    if latest_ago is not None and latest_ago >= NO_ACTIVITY_WARNING_SECONDS:
        print(
            f"  warning: no observed output file modification in the last {int(latest_ago // 60)} minutes",
            flush=True,
        )
        print("  process is still alive; check logs or system monitor if concerned", flush=True)
    print("  full logs:", flush=True)
    print(f"    stdout: {stdout_log}", flush=True)
    print(f"    stderr: {stderr_log}", flush=True)
    if tail_log_lines > 0:
        for label, log in (("stdout", stdout_log), ("stderr", stderr_log)):
            lines = useful_tail_lines(log, n=tail_log_lines, max_chars=max_progress_log_line_chars)
            if lines:
                print(f"  recent {label}:", flush=True)
                for line in lines:
                    print(f"    {line}", flush=True)


def _pump_stream(stream, log_fh) -> None:
    try:
        for line in iter(stream.readline, ""):
            if line == "":
                break
            log_fh.write(line)
            log_fh.flush()
    finally:
        try:
            stream.close()
        except Exception:
            pass


def run_logged_command_with_progress(
    command: Sequence[str],
    *,
    cwd: Path,
    stdout_log: Path,
    stderr_log: Path,
    progress_label: str,
    progress_interval_seconds: int = 60,
    quiet: bool = False,
    watch_paths: Sequence[Path] = (),
    tail_log_lines: int = 0,
    max_progress_log_line_chars: int = 300,
    progress_mode: str = "normal",
) -> int:
    """Run a command with compact production-friendly progress output.

    Full stdout/stderr are always preserved in log files. Terminal log tails are
    disabled by default and, when explicitly enabled, are filtered/truncated.
    """
    stdout_log.parent.mkdir(parents=True, exist_ok=True)
    stderr_log.parent.mkdir(parents=True, exist_ok=True)
    progress_interval_seconds = max(1, int(progress_interval_seconds or 60))
    tail_log_lines = max(0, int(tail_log_lines or 0))
    max_progress_log_line_chars = max(40, int(max_progress_log_line_chars or 300))
    if progress_mode not in {"compact", "normal", "verbose"}:
        progress_mode = "normal"

    with stdout_log.open("w", encoding="utf-8") as out_fh, stderr_log.open("w", encoding="utf-8") as err_fh:
        proc = subprocess.Popen(
            list(map(str, command)),
            cwd=str(cwd),
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            bufsize=1,
        )
        assert proc.stdout is not None
        assert proc.stderr is not None
        stdout_thread = threading.Thread(target=_pump_stream, args=(proc.stdout, out_fh), daemon=True)
        stderr_thread = threading.Thread(target=_pump_stream, args=(proc.stderr, err_fh), daemon=True)
        stdout_thread.start()
        stderr_thread.start()

        start_time = time.time()
        last_report = start_time
        heartbeat_count = 0
        if not quiet:
            print(f"[{progress_label}] started", flush=True)
            print(f"  command: {' '.join(map(str, command))}", flush=True)
            print(f"  cwd: {cwd}", flush=True)
            print("  logs:", flush=True)
            print(f"    stdout: {stdout_log}", flush=True)
            print(f"    stderr: {stderr_log}", flush=True)
            print("", flush=True)

        while proc.poll() is None:
            now = time.time()
            if not quiet and now - last_report >= progress_interval_seconds:
                heartbeat_count += 1
                if progress_mode == "compact":
                    print_compact_progress(
                        progress_label=progress_label,
                        stderr_log=stderr_log,
                        start_time=start_time,
                        watch_paths=watch_paths,
                    )
                elif progress_mode == "verbose":
                    print_expanded_progress(
                        progress_label=progress_label,
                        stderr_log=stderr_log,
                        stdout_log=stdout_log,
                        start_time=start_time,
                        watch_paths=watch_paths,
                        tail_log_lines=tail_log_lines,
                        max_progress_log_line_chars=max_progress_log_line_chars,
                        show_paths=True,
                        command=command,
                        cwd=cwd,
                    )
                else:
                    print_compact_progress(
                        progress_label=progress_label,
                        stderr_log=stderr_log,
                        start_time=start_time,
                        watch_paths=watch_paths,
                    )
                    # Every fifth heartbeat, include a compact expanded status block.
                    if heartbeat_count % 5 == 0:
                        print_expanded_progress(
                            progress_label=progress_label,
                            stderr_log=stderr_log,
                            stdout_log=stdout_log,
                            start_time=start_time,
                            watch_paths=watch_paths,
                            tail_log_lines=tail_log_lines,
                            max_progress_log_line_chars=max_progress_log_line_chars,
                            show_paths=False,
                            command=command,
                            cwd=cwd,
                        )
                last_report = now
            time.sleep(2)

        stdout_thread.join(timeout=5)
        stderr_thread.join(timeout=5)

        if not quiet:
            print(
                f"[{progress_label}] finished after {seconds_to_hhmmss(time.time() - start_time)} with return code {proc.returncode}",
                flush=True,
            )

        return int(proc.returncode or 0)


def observed_sizes(paths_by_name: dict[str, Path]) -> dict[str, int]:
    return {name: directory_size_bytes(path) for name, path in paths_by_name.items()}


def observed_mtimes(paths_by_name: dict[str, Path]) -> dict[str, Optional[str]]:
    return {name: mtime_to_utc(latest_mtime(path)) for name, path in paths_by_name.items()}


def write_manifest(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
