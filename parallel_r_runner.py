#!/usr/bin/env python3
"""
Parallel executor for R code using rpy2 with multiprocessing + threading.

Key features:
- Spawn-safe R initialization per worker (uses spawn, not fork).
- Batch-level isolation with per-task environment injection.
- Live progress monitor with ETA, throughput, elapsed time.
- Resilient error reporting without blocking other tasks.

Prerequisites:
- Python 3.9+.
- rpy2 installed and able to load your R installation.
- The R script should expose an entrypoint function that can be called
  programmatically (see the __main__ block for an example wrapper).
"""

from __future__ import annotations

import argparse
import os
import sys
import time
import traceback
from collections import deque
from dataclasses import dataclass, field
from multiprocessing import get_context
from queue import Empty
from threading import Event, Thread
from typing import Any, Dict, Iterable, List, Optional, Sequence


# Default base directory for example batches (mirrors the R config snippet)
DEFAULT_EXAMPLE_BASE_DIR = "/home/Dementia_Depression_MR-analysis/data/MR_pipeline_demo"

# Track whether the R script has been sourced inside a worker process.
# This avoids re-sourcing (and re-initializing caches) for every batch.
_R_SESSION_STATE = {"script_path": None}

# Global progress queue injected into worker processes via initializer.
_PROGRESS_QUEUE = None

# Cache for completion markers recovered from parent-level logs / manifests.
_COMPLETED_BATCH_RECORDS: Dict[str, Dict[str, set[str]]] = {}


@dataclass
class RRuntimeConfig:
    """Configuration for R execution inside each worker."""

    script_path: str
    entrypoint: Optional[str] = "run_main_analysis"  # R function to call after sourcing
    r_options: Dict[str, Any] = field(default_factory=dict)  # Passed to options()
    env: Dict[str, str] = field(default_factory=dict)  # Extra env vars for R
    quiet: bool = True  # Suppress R console output if True
    dt_threads: Optional[int] = None  # data.table thread count inside R (setDTthreads)
    r_workers: Optional[int] = None  # Controls mr_parallel_workers inside R


@dataclass
class BatchSpec:
    """
    Description of a single batch task.
    parameters: passed as keyword arguments to the R entrypoint.
    """

    name: str
    parameters: Dict[str, Any]
    output_dir: str
    env: Dict[str, str] = field(default_factory=dict)  # Per-batch env overrides


def _ensure_output_dir(path: str) -> None:
    os.makedirs(path, exist_ok=True)


def _env_bool(name: str, default: bool = False) -> bool:
    val = os.environ.get(name)
    if val is None:
        return default
    return val.strip().lower() in {"1", "true", "yes", "y", "on"}


def _env_int(name: str) -> Optional[int]:
    val = os.environ.get(name)
    if val is None or val == "":
        return None
    try:
        return int(val)
    except ValueError:
        return None


def _file_size(path: str) -> Optional[int]:
    try:
        return os.stat(path).st_size
    except OSError:
        return None


def _effective_processes(requested: Optional[int], batch_count: int) -> int:
    cpu_count = os.cpu_count() or 1
    if requested is None or requested <= 0:
        requested = cpu_count
    return max(1, min(requested, max(1, batch_count)))


def _extract_logged_result_path(line: str) -> Optional[str]:
    """
    Extract the logged result CSV path from a task log line.

    Expected format:
        2026-01-01 12:34:56<TAB>uvmr_results.csv /path/to/file.csv rows 5
    """
    if "\t" not in line:
        return None

    _, message = line.rstrip("\n").split("\t", 1)
    parts = message.strip().split()
    if len(parts) < 2:
        return None

    result_name = parts[0].lower()
    if result_name not in {
        "uvmr_results.csv",
        "mvmr_results.csv",
        "mediation_results.csv",
    }:
        return None

    return parts[1]


def _load_completed_batch_records(parent_dir: str) -> Dict[str, set[str]]:
    """
    Recover completed batch markers from lightweight parent-level metadata.

    This supports compacting finished batch directories after their outputs have
    been merged elsewhere. The runner can still skip them later via:
    1. completed_batches.txt
    2. task logs such as uvmr.log / mvmr.log / mediation.log
    """
    normalized_parent = os.path.normpath(parent_dir)
    cached = _COMPLETED_BATCH_RECORDS.get(normalized_parent)
    if cached is not None:
        return cached

    batch_names: set[str] = set()
    batch_dirs: set[str] = set()

    manifest_path = os.path.join(normalized_parent, "completed_batches.txt")
    if os.path.exists(manifest_path):
        try:
            with open(manifest_path, "r", encoding="utf-8") as handle:
                for raw_line in handle:
                    line = raw_line.strip()
                    if not line or line.startswith("#"):
                        continue

                    entry = line.split("\t", 1)[0].strip()
                    if not entry:
                        continue

                    batch_name = os.path.basename(os.path.normpath(entry))
                    if batch_name:
                        batch_names.add(batch_name)
                    if os.path.isabs(entry):
                        batch_dirs.add(os.path.normpath(entry))
        except (OSError, UnicodeDecodeError):
            pass

    for log_name in ("uvmr.log", "mvmr.log", "mediation.log"):
        log_path = os.path.join(normalized_parent, log_name)
        if not os.path.exists(log_path):
            continue

        try:
            with open(log_path, "r", encoding="utf-8") as handle:
                for raw_line in handle:
                    csv_path = _extract_logged_result_path(raw_line)
                    if not csv_path:
                        continue

                    batch_dir = os.path.normpath(os.path.dirname(csv_path))
                    batch_name = os.path.basename(batch_dir)
                    if batch_name:
                        batch_names.add(batch_name)
                    if batch_dir:
                        batch_dirs.add(batch_dir)
        except (OSError, UnicodeDecodeError):
            continue

    records = {"names": batch_names, "dirs": batch_dirs}
    _COMPLETED_BATCH_RECORDS[normalized_parent] = records
    return records


def _is_batch_completed_via_parent_records(output_dir: str) -> bool:
    """
    Check whether a batch is recorded as completed in parent-level metadata.
    """
    normalized_output_dir = os.path.normpath(output_dir)
    parent_dir = os.path.dirname(normalized_output_dir)
    batch_name = os.path.basename(normalized_output_dir)
    records = _load_completed_batch_records(parent_dir)
    return (
        normalized_output_dir in records["dirs"]
        or batch_name in records["names"]
    )


def _ensure_r_sourced(ro, config: "RRuntimeConfig") -> None:
    """
    Source the R script only once per worker process (per script path).
    This preserves in-session caches and reduces overhead.
    """
    global _R_SESSION_STATE
    needs_source = _R_SESSION_STATE.get("script_path") != config.script_path
    if not needs_source and config.entrypoint:
        try:
            needs_source = config.entrypoint not in ro.globalenv
        except Exception:  # noqa: BLE001
            needs_source = True
    if needs_source:
        ro.r["source"](config.script_path)
        _R_SESSION_STATE["script_path"] = config.script_path


def _init_worker(progress_queue) -> None:
    """Initializer for worker processes to receive the progress queue."""
    global _PROGRESS_QUEUE
    _PROGRESS_QUEUE = progress_queue


def _is_batch_completed(output_dir: str) -> bool:
    """
    检查批次是否已经完成。

    判断标准（按优先级）：
    1. progress.log 中存在 "run_end" 标记（最可靠）
    2. 存在 .completed 标记文件
    3. 存在主要的结果文件（UVMR_results.csv, MVMR_results.csv, Mediation_results.csv）
    4. 父目录中的 completed_batches.txt 或任务日志记录了该批次已完成

    Returns:
        True if batch is completed, False otherwise
    """
    recorded_completed = _is_batch_completed_via_parent_records(output_dir)

    if not os.path.isdir(output_dir):
        return recorded_completed

    # 优先检查 progress.log 中的 run_end 标记（最可靠的完成标志）
    progress_log = os.path.join(output_dir, "progress.log")
    if os.path.exists(progress_log):
        try:
            with open(progress_log, "r", encoding="utf-8") as f:
                for line in f:
                    # 格式: PROGRESS\trun_end\t1\t1\t时间戳
                    if line.startswith("PROGRESS") and "run_end" in line:
                        parts = line.strip().split("\t")
                        if len(parts) >= 2 and parts[1] == "run_end":
                            return True
        except (OSError, UnicodeDecodeError):
            pass  # 如果读取失败，继续检查其他标记

    # 检查完成标记文件
    completed_marker = os.path.join(output_dir, ".completed")
    if os.path.exists(completed_marker):
        return True

    # 检查是否存在主要结果文件（至少有一个）
    result_files = [
        "UVMR_results.csv",
        "MVMR_results.csv",
        "Mediation_results.csv",
        "uvmr_results.csv",
        "mvmr_results.csv",
        "mediation_results.csv"
    ]

    for result_file in result_files:
        result_path = os.path.join(output_dir, result_file)
        if os.path.exists(result_path) and os.path.getsize(result_path) > 0:
            return True

    return recorded_completed


def _mark_batch_completed(output_dir: str) -> None:
    """
    标记批次为已完成。
    在输出目录创建 .completed 文件。
    """
    try:
        _ensure_output_dir(output_dir)
        completed_marker = os.path.join(output_dir, ".completed")
        with open(completed_marker, "w", encoding="utf-8") as f:
            f.write(f"Completed at: {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
    except OSError:
        pass  # 如果无法创建标记文件，不影响主流程


def _prepare_tmpdir(base: str, batch_name: str) -> str:
    """
    Create an isolated temp directory per batch to avoid /tmp exhaustion
    or permissions issues (addresses 'Fatal error: creating temporary file for -e failed').
    """
    tmp_dir = os.path.join(base, "tmp", batch_name)
    os.makedirs(tmp_dir, exist_ok=True)
    return tmp_dir


class _silence_r_console:
    """
    Temporarily silence R console output via rpy2 callbacks.
    """

    def __init__(self, enabled: bool = True):
        self.enabled = enabled
        self._orig_print = None
        self._orig_warn = None
        self._warn_attr = None

    def __enter__(self):
        if not self.enabled:
            return
        from rpy2.rinterface_lib import callbacks

        self._orig_print = callbacks.consolewrite_print
        # rpy2 versions differ: consolewrite_warn or consolewrite_warnerror
        if hasattr(callbacks, "consolewrite_warn"):
            self._warn_attr = "consolewrite_warn"
        elif hasattr(callbacks, "consolewrite_warnerror"):
            self._warn_attr = "consolewrite_warnerror"
        else:
            self._warn_attr = None

        callbacks.consolewrite_print = lambda x: None
        if self._warn_attr:
            self._orig_warn = getattr(callbacks, self._warn_attr)
            setattr(callbacks, self._warn_attr, lambda x: None)

    def __exit__(self, exc_type, exc, tb):
        if not self.enabled:
            return
        from rpy2.rinterface_lib import callbacks

        callbacks.consolewrite_print = self._orig_print
        if self._warn_attr and self._orig_warn is not None:
            setattr(callbacks, self._warn_attr, self._orig_warn)


def _list_gwas_files(dir_path: str) -> List[str]:
    """
    Mirror of the R helper in Demo_Test_Analysis.R: list GWAS files in a directory.
    """
    if not os.path.isdir(dir_path):
        return []
    return sorted(
        [
            os.path.join(dir_path, f)
            for f in os.listdir(dir_path)
            if any(
                f.lower().endswith(ext)
                for ext in (
                    ".tsv",
                    ".tsv.gz",
                    ".gz.tsv",
                    ".txt.gz",
                    ".txt",
                    ".gz.txt",
                )
            )
        ]
    )


def _select_demo_files(base_dir: str) -> Dict[str, List[str]]:
    """
    Build a demo subset following Demo_Test_Analysis.R:
      - 1 exposure from each exposure folder (up to 3 total)
      - 2 mediators
      - 2 outcomes
      - all covariates
    """
    expo_dirs = [
        os.path.join(base_dir, "Standardized Circulating human plasma proteome_Data"),
        os.path.join(base_dir, "Standardized Circulating metabolic biomarkers_Data"),
        os.path.join(base_dir, "Circulating inflammatory proteins_Data"),
    ]
    medi_dir = os.path.join(base_dir, "Cerebrospinal fluid metabolomics_Data")
    out_dir = os.path.join(base_dir, "Outcomes")
    covar_dir = os.path.join(base_dir, "Covariates_SES")

    exposures: List[str] = []
    for d in expo_dirs:
        files = _list_gwas_files(d)
        if files:
            exposures.append(files[0])

    mediators = _list_gwas_files(medi_dir)[:2]
    outcomes = _list_gwas_files(out_dir)[:2]
    covariates = _list_gwas_files(covar_dir)

    return {
        "exposures": exposures,
        "mediators": mediators,
        "outcomes": outcomes,
        "covariates": covariates,
    }


def _execute_r(batch: BatchSpec, config: RRuntimeConfig) -> Dict[str, Any]:
    """
    Run one batch in R via rpy2.

    Notes:
    - rpy2 is imported inside the worker to avoid fork-after-init issues.
    - Uses spawn context; do not import rpy2 at module import time in the parent.
    - Returns metadata only (no R objects) to avoid pickle/segfault issues.
    """
    # Late import inside worker process
    import rpy2.robjects as ro
    from rpy2.robjects import pandas2ri

    # Apply environment overrides for this batch
    env = dict(os.environ)
    env.update(config.env)
    env.update(batch.env)
    os.environ.update(env)

    # Configure R session
    if config.r_options:
        ro.r.options(**config.r_options)

    with _silence_r_console(enabled=config.quiet):
        # Optionally force R worker count via option mr_parallel_workers
        if config.r_workers:
            ro.r(f"options(mr_parallel_workers = {int(config.r_workers)})")

        # Tune data.table threading if requested
        if config.dt_threads:
            ro.r(f"data.table::setDTthreads({int(config.dt_threads)})")

        # Expose batch name for logging on the R side if desired
        ro.globalenv[".py_batch_name"] = batch.name

        # Source the R script once per worker process
        _ensure_r_sourced(ro, config)

        # Call entrypoint if provided
        result_meta = {"entrypoint_called": False, "output_files": []}

        if config.entrypoint:
            if config.entrypoint not in ro.globalenv:
                raise RuntimeError(
                    f"Entrypoint '{config.entrypoint}' not found in R after sourcing {config.script_path}"
                )
            r_fun = ro.globalenv[config.entrypoint]

            # Convert Python types (including pandas) to R
            with ro.conversion.localconverter(ro.default_converter + pandas2ri.converter):
                # Call R function but don't return the R object (causes segfault)
                r_fun(**batch.parameters, output_dir=batch.output_dir)

            result_meta["entrypoint_called"] = True

            # Collect output files (safe to return)
            if os.path.isdir(batch.output_dir):
                result_meta["output_files"] = [
                    f for f in os.listdir(batch.output_dir)
                    if f.endswith(('.csv', '.tsv', '.txt', '.log'))
                ]

        # Encourage R to free memory within the worker
        ro.r("gc(verbose = FALSE)")

        return result_meta


def _worker(
    batch: BatchSpec,
    config: RRuntimeConfig,
    progress_queue=None,
    progress_log: bool = True,
) -> Dict[str, Any]:
    """
    Worker wrapper to execute a batch and post progress events.
    Includes robust error handling to prevent segfaults from crashing the pool.
    """
    if progress_queue is None:
        progress_queue = _PROGRESS_QUEUE

    def _emit(msg: Dict[str, Any]) -> None:
        if progress_queue is None:
            return
        try:
            progress_queue.put(msg)
        except Exception:  # noqa: BLE001
            pass
    start = time.monotonic()
    _emit({"event": "start", "batch": batch.name})

    # Setup isolation and error state
    result_dict = {
        "batch": batch.name,
        "status": "error",
        "error": "Unknown error",
        "traceback": "",
        "duration": 0.0,
    }

    try:
        _ensure_output_dir(batch.output_dir)
        progress_file = os.path.join(batch.output_dir, "progress.log")
        error_file = os.path.join(batch.output_dir, "error.log")
        run_log_file = os.path.join(batch.output_dir, "run.log")
        try:
            if progress_log:
                with open(progress_file, "w", encoding="utf-8"):
                    pass  # truncate if exists
            with open(error_file, "w", encoding="utf-8"):
                pass  # truncate if exists
            with open(run_log_file, "w", encoding="utf-8"):
                pass  # truncate if exists
        except OSError:
            pass

        # Inject per-batch progress, error, and run log files for R to write into
        batch.env = dict(batch.env)
        if progress_log:
            batch.env["PY_PROGRESS_FILE"] = progress_file
        batch.env["PY_ERROR_LOG"] = error_file
        batch.env["PY_RUN_LOG"] = run_log_file
        tmp_dir = _prepare_tmpdir(batch.output_dir, batch.name)

        # Isolate temp directories to avoid permission issues (/tmp) and contention
        os.environ["TMPDIR"] = tmp_dir
        os.environ["TEMP"] = tmp_dir
        os.environ["TEMPDIR"] = tmp_dir
        os.environ["R_TMPDIR"] = tmp_dir

        # Execute R code
        result_meta = _execute_r(batch, config)

        # Success path
        duration = time.monotonic() - start

        # 标记批次为已完成
        _mark_batch_completed(batch.output_dir)

        _emit(
            {"event": "done", "batch": batch.name, "status": "ok", "duration": duration}
        )
        result_dict = {
            "batch": batch.name,
            "status": "ok",
            "result": result_meta,
            "duration": duration,
        }

    except SystemExit:
        # Catch SystemExit to prevent worker from dying silently
        duration = time.monotonic() - start
        error_msg = "R script called quit() or stop()"
        _emit(
            {
                "event": "done",
                "batch": batch.name,
                "status": "error",
                "duration": duration,
                "error": error_msg,
            }
        )
        result_dict["error"] = error_msg
        result_dict["duration"] = duration
        result_dict["traceback"] = "R SystemExit"

    except MemoryError:
        # Out of memory
        duration = time.monotonic() - start
        error_msg = "Out of memory (R or Python)"
        _emit(
            {
                "event": "done",
                "batch": batch.name,
                "status": "error",
                "duration": duration,
                "error": error_msg,
            }
        )
        result_dict["error"] = error_msg
        result_dict["duration"] = duration
        result_dict["traceback"] = traceback.format_exc()

    except Exception as exc:  # noqa: BLE001
        # General exceptions
        duration = time.monotonic() - start
        _emit(
            {
                "event": "done",
                "batch": batch.name,
                "status": "error",
                "duration": duration,
                "error": str(exc),
            }
        )
        result_dict["error"] = str(exc)
        result_dict["duration"] = duration
        result_dict["traceback"] = traceback.format_exc()

    finally:
        # Force garbage collection to free R resources
        try:
            import gc
            gc.collect()
        except:  # noqa: E722
            pass

    return result_dict


def _progress_monitor(
    queue,
    stop_event: Event,
    batch_names,
    progress_files,
    error_files,
    refresh: float = 0.5,
    io_refresh: Optional[float] = None,
) -> None:
    """
    Live progress monitor running in a dedicated thread.
    Clears the screen and redraws a per-batch table with ETA.
    Also monitors error.log files for MRlap and other errors.
    io_refresh controls how often progress/error logs are polled.
    """
    total_batches = max(1, len(batch_names))
    started = 0
    finished = 0
    failed = 0
    in_flight = 0
    start_time = time.monotonic()
    io_refresh = refresh if io_refresh is None else max(0.1, float(io_refresh))
    last_poll = 0.0
    statuses = {name: {"status": "pending", "duration": None} for name in batch_names}
    batch_order = {name: i for i, name in enumerate(batch_names)}
    running_set = set()
    running_progress = {}
    recent_completed = deque(maxlen=5)
    failed_batches = []
    start_times = {}
    progress = {name: "" for name in batch_names}
    offsets = {name: 0 for name in batch_names}
    error_offsets = {name: 0 for name in batch_names}
    progress_sizes = {}
    error_sizes = {}
    last_logged = {name: "" for name in batch_names}
    errors_logged = {name: [] for name in batch_names}

    def parse_fraction(text: str) -> float:
        """
        Parse the last a/b fragment in text and return a/b as float in [0,1].
        Handles edge cases where a > b (returns 1.0) or malformed fractions.
        """
        if not text:
            return 0.0
        parts = text.strip().split()
        for token in reversed(parts):
            if "/" in token:
                nums = token.split("/", 1)
                try:
                    a = float(nums[0])
                    b = float(nums[1])
                    if b <= 0:
                        return 0.0
                    # Clamp to [0, 1] - if a > b, treat as complete (1.0)
                    ratio = a / b
                    return max(0.0, min(1.0, ratio))
                except (ValueError, ZeroDivisionError):
                    continue  # Try next token instead of returning 0
        return 0.0

    def calculate_batch_progress(progress_text: str) -> float:
        """
        Calculate weighted progress for a batch based on the current step.
        R script has multiple sequential steps, each weighted by importance.

        Steps and weights:
        - uvmr_clump/uvmr_load: 0.0 - 0.05 (5% - data loading)
        - uvmr_pairs: 0.05 - 0.25 (20% - UVMR analysis)
        - mvmr_clump/mvmr_load: 0.25 - 0.30 (5% - data loading)
        - mvmr_start/mvmr_pairs: 0.30 - 0.50 (20% - MVMR analysis)
        - mediation_clump/mediation_load: 0.50 - 0.55 (5% - data loading)
        - mediation_start/mediation_triplets: 0.55 - 1.0 (45% - Mediation analysis)
        """
        if not progress_text:
            return 0.0

        parts = progress_text.strip().split()
        if not parts:
            return 0.0

        step = parts[0]
        frac = parse_fraction(progress_text)

        # Map step to weighted progress range
        if step.startswith("uvmr_clump") or step.startswith("uvmr_load"):
            return 0.0 + (frac * 0.05)  # 0-5%
        elif step.startswith("uvmr"):
            return 0.05 + (frac * 0.20)  # 5-25%
        elif step.startswith("mvmr_clump") or step.startswith("mvmr_load"):
            return 0.25 + (frac * 0.05)  # 25-30%
        elif step.startswith("mvmr"):
            return 0.30 + (frac * 0.20)  # 30-50%
        elif step.startswith("mediation_clump") or step.startswith("mediation_load"):
            return 0.50 + (frac * 0.05)  # 50-55%
        elif step.startswith("mediation"):
            return 0.55 + (frac * 0.45)  # 55-100%
        elif step.startswith("run_"):
            # Special markers for overall run status
            return frac  # Use as-is
        else:
            # Unknown step, use raw fraction
            return frac

    def format_progress_text(step: str, done: Optional[int] = None, total: Optional[int] = None, msg: str = "") -> str:
        text = step
        if done is not None and total is not None:
            text = f"{step} {done}/{total}"
        msg = (msg or "").strip().replace("\n", " ")
        if msg:
            if len(msg) > 40:
                msg = msg[:37] + "..."
            text = f"{text} | {msg}"
        return text

    def fmt_eta(rate: float, remaining: int) -> str:
        if rate <= 0:
            return "--:--:--"
        seconds = remaining / rate
        return time.strftime("%H:%M:%S", time.gmtime(seconds))

    def render():
        elapsed = max(1e-6, time.monotonic() - start_time)
        # Compute effective progress including running batches with partial fractions
        # Use weighted progress calculation to account for multi-step R processing
        running_sum = sum(running_progress.values()) if running_progress else 0.0
        eff_progress = (finished + running_sum) / max(1, total_batches)
        pct = eff_progress * 100
        rate = eff_progress / elapsed if elapsed > 0 else 0.0
        eta = fmt_eta(rate, 1 - eff_progress) if eff_progress > 0 else "--:--:--"

        sys.stdout.write("\x1b[2J\x1b[H")  # clear screen
        sys.stdout.write(
            f"Batches: {finished}/{total_batches} ({pct:5.1f}%) | "
            f"running={in_flight} | failed={failed} | "
            f"elapsed={elapsed:6.1f}s | eta={eta}\n"
        )
        sys.stdout.write("=" * 80 + "\n")

        # For large batch counts, only show running/recent batches to avoid screen overflow
        if total_batches > 20:
            # Show only running batches and last few completed/failed
            running_batches = sorted(running_set, key=batch_order.get)

            sys.stdout.write("Running batches:\n")
            sys.stdout.write("-" * 80 + "\n")
            if running_batches:
                for name in running_batches:
                    info = statuses.get(name, {})
                    run_for = time.monotonic() - start_times.get(name, start_time)
                    dur_str = f"{run_for:6.1f}s"
                    prog = progress.get(name, "")
                    sys.stdout.write(f"{name:25s} | running | {dur_str} | {prog}\n")
            else:
                sys.stdout.write("  (none)\n")

            # Show last 5 completed
            if recent_completed:
                sys.stdout.write("\nRecent completed (last 5):\n")
                sys.stdout.write("-" * 80 + "\n")
                for name, info in list(recent_completed):
                    dur = info.get("duration", 0)
                    dur_str = f"{dur:6.1f}s"
                    sys.stdout.write(f"{name:25s} | ok      | {dur_str}\n")

            # Show all failed
            if failed_batches:
                sys.stdout.write("\nFailed batches:\n")
                sys.stdout.write("-" * 80 + "\n")
                for name, info in failed_batches:
                    dur = info.get("duration", 0)
                    dur_str = f"{dur:6.1f}s"
                    sys.stdout.write(f"{name:25s} | error   | {dur_str}\n")
        else:
            # Small batch count: show all batches as before
            sys.stdout.write("All batches:\n")
            sys.stdout.write("-" * 80 + "\n")
            for name in batch_names:
                info = statuses.get(name, {})
                st = info.get("status", "pending")
                dur = info.get("duration")
                if st == "running":
                    run_for = time.monotonic() - start_times.get(name, start_time)
                    dur_str = f"{run_for:6.1f}s"
                elif dur is not None:
                    dur_str = f"{dur:6.1f}s"
                else:
                    dur_str = "   --  "
                prog = progress.get(name, "")
                sys.stdout.write(f"{name:25s} | {st:7s} | {dur_str} | {prog}\n")

        sys.stdout.flush()

    last_render = 0.0
    while not stop_event.is_set() or not queue.empty():
        render_needed = False
        try:
            msg = queue.get(timeout=refresh)
        except Empty:
            msg = None

        if msg:
            if msg["event"] == "start":
                started += 1
                in_flight += 1
                statuses[msg["batch"]] = {"status": "running", "duration": None}
                running_set.add(msg["batch"])
                running_progress[msg["batch"]] = 0.0
                start_times[msg["batch"]] = time.monotonic()
                sys.stdout.write(f"[start] {msg['batch']} ({started}/{total_batches})\n")
                sys.stdout.flush()
                render_needed = True
            elif msg["event"] == "done":
                finished += 1
                in_flight = max(0, in_flight - 1)
                if msg.get("status") == "error":
                    failed += 1
                statuses[msg["batch"]] = {
                    "status": msg.get("status", "done"),
                    "duration": msg.get("duration", 0.0),
                }
                running_set.discard(msg["batch"])
                running_progress.pop(msg["batch"], None)
                if msg.get("status") == "error":
                    failed_batches.append((msg["batch"], statuses[msg["batch"]]))
                else:
                    recent_completed.append((msg["batch"], statuses[msg["batch"]]))
                sys.stdout.write(
                    f"[done] {msg['batch']} {msg.get('status', 'done')} "
                    f"({msg.get('duration', 0.0):.1f}s)\n"
                )
                sys.stdout.flush()
                render_needed = True

        now = time.monotonic()
        if (now - last_poll) >= io_refresh:
            last_poll = now
            # Only poll logs for running batches to reduce I/O
            running_batches = list(running_set)

            # Read incremental progress updates written by R
            for name in running_batches:
                path = progress_files.get(name)
                if not path:
                    continue
                size = _file_size(path)
                if size is None:
                    continue
                if progress_sizes.get(name) == size:
                    continue
                progress_sizes[name] = size
                if size < offsets.get(name, 0):
                    offsets[name] = 0
                try:
                    with open(path, "r", encoding="utf-8") as fh:
                        fh.seek(offsets.get(name, 0))
                        lines = fh.readlines()
                        offsets[name] = fh.tell()
                    for line in lines:
                        parts = line.strip().split("\t")
                        if len(parts) >= 4 and parts[0] == "PROGRESS":
                            step = parts[1]
                            msg = parts[4] if len(parts) > 4 else ""
                            try:
                                done = int(parts[2])
                                total = int(parts[3])
                                progress[name] = format_progress_text(step, done, total, msg)
                            except ValueError:
                                progress[name] = format_progress_text(step, msg=msg)
                            running_progress[name] = calculate_batch_progress(progress[name])
                            if progress[name] != last_logged.get(name, ""):
                                sys.stdout.write(f"[progress] {name} {progress[name]}\n")
                                sys.stdout.flush()
                                last_logged[name] = progress[name]
                                render_needed = True
                except OSError:
                    continue

            # Read incremental error logs written by R
            for name in running_batches:
                path = error_files.get(name)
                if not path:
                    continue
                size = _file_size(path)
                if size is None:
                    continue
                if error_sizes.get(name) == size:
                    continue
                error_sizes[name] = size
                if size < error_offsets.get(name, 0):
                    error_offsets[name] = 0
                try:
                    with open(path, "r", encoding="utf-8") as fh:
                        fh.seek(error_offsets.get(name, 0))
                        lines = fh.readlines()
                        error_offsets[name] = fh.tell()
                    for line in lines:
                        parts = line.strip().split("\t")
                        # Format: ERROR | timestamp | error_type | exposure | outcome | location | details
                        if len(parts) >= 3 and parts[0] == "ERROR":
                            error_type = parts[2] if len(parts) > 2 else "Unknown"
                            details = parts[6] if len(parts) > 6 else ""
                            error_msg = f"[error] {name} {error_type}: {details}"
                            if error_msg not in errors_logged[name]:
                                sys.stdout.write(error_msg + "\n")
                                sys.stdout.flush()
                                errors_logged[name].append(error_msg)
                                render_needed = True
                except OSError:
                    continue
        if render_needed or (now - last_render) >= refresh:
            render()
            last_render = now

    # Final flush of progress updates after stop_event is set
    for name, path in progress_files.items():
        try:
            if not os.path.exists(path):
                continue
            with open(path, "r", encoding="utf-8") as fh:
                fh.seek(offsets.get(name, 0))
                lines = fh.readlines()
            for line in lines:
                parts = line.strip().split("\t")
                if len(parts) >= 4 and parts[0] == "PROGRESS":
                    step = parts[1]
                    msg = parts[4] if len(parts) > 4 else ""
                    try:
                        done = int(parts[2])
                        total = int(parts[3])
                        progress[name] = format_progress_text(step, done, total, msg)
                    except ValueError:
                        progress[name] = format_progress_text(step, msg=msg)
                    if progress[name] != last_logged.get(name, ""):
                        sys.stdout.write(f"[progress] {name} {progress[name]}\n")
                        sys.stdout.flush()
                        last_logged[name] = progress[name]
        except OSError:
            continue

    render()
    sys.stdout.write("\n")
    sys.stdout.flush()


def run_batches(
    batches: Iterable[BatchSpec],
    config: RRuntimeConfig,
    processes: Optional[int] = None,
    skip_completed: bool = False,
    monitor: bool = True,
    monitor_refresh: float = 0.5,
    monitor_io: Optional[float] = None,
    progress_log: bool = True,
) -> List[Dict[str, Any]]:
    """
    Execute batches in parallel using multiprocessing and a threaded progress monitor.

    Args:
        batches: Iterable of BatchSpec objects to process
        config: R runtime configuration
        processes: Number of worker processes (defaults to CPU count)
        skip_completed: If True, skip batches that have already been completed
        monitor: If True, show live progress monitor
        monitor_refresh: Refresh interval for the progress monitor (seconds)
        monitor_io: Polling interval for progress/error logs (seconds)
        progress_log: Whether to enable progress.log writing by R
    """
    batch_list = list(batches)

    # 过滤已完成的批次
    if skip_completed:
        original_count = len(batch_list)
        skipped_batches = []
        pending_batches = []

        for batch in batch_list:
            if _is_batch_completed(batch.output_dir):
                skipped_batches.append(batch.name)
            else:
                pending_batches.append(batch)

        batch_list = pending_batches

        if skipped_batches:
            print(f"\n[SKIP] Skipping {len(skipped_batches)} completed batches (out of {original_count} total)")
            print(f"[SKIP] Processing {len(batch_list)} remaining batches")
            if len(skipped_batches) <= 10:
                print("[SKIP] Skipped batches:")
                for name in skipped_batches:
                    print(f"  - {name}")
            else:
                print(f"[SKIP] First 5 skipped batches:")
                for name in skipped_batches[:5]:
                    print(f"  - {name}")
                print(f"  ... and {len(skipped_batches) - 5} more")
            print()

        if not batch_list:
            print("[INFO] All batches have been completed. Nothing to process.")
            return []
    # Respect requested worker count; default to CPU count and never exceed batch count
    if processes is None or processes <= 0:
        processes = os.cpu_count() or 1
    processes = max(1, min(processes, len(batch_list)))
    ctx = get_context("spawn")  # R is not fork-safe; spawn is safest
    progress_queue = ctx.Queue() if monitor else None
    stop_event = Event()
    if monitor:
        progress_files = (
            {b.name: os.path.join(b.output_dir, "progress.log") for b in batch_list}
            if progress_log else {}
        )
        error_files = {b.name: os.path.join(b.output_dir, "error.log") for b in batch_list}

        monitor_thread = Thread(
            target=_progress_monitor,
            args=(
                progress_queue,
                stop_event,
                [b.name for b in batch_list],
                progress_files,
                error_files,
                monitor_refresh,
                monitor_io,
            ),
            daemon=True,
        )
        monitor_thread.start()

    results: List[Any] = []
    with ctx.Pool(processes=processes, initializer=_init_worker, initargs=(progress_queue,)) as pool:
        async_results = [
            pool.apply_async(_worker, (batch, config, None, progress_log))
            for batch in batch_list
        ]
        pool.close()
        pool.join()
        results = [res.get() for res in async_results]

    if monitor:
        stop_event.set()
        monitor_thread.join()
    return results


def _build_example_batches(base_dir: str) -> List[BatchSpec]:
    """
    Example helper showing how to define batches.

    Follows the selection logic in Demo_Test_Analysis.R:
      - 1 exposure from each exposure folder (up to 3 total)
      - 2 mediators (round-robin across batches)
      - 2 outcomes
      - all covariates passed as a list (if the R entrypoint supports it)

    Assumes the R entrypoint signature:
        run_main_analysis(exposure_path, mediator_path, outcome_path,
                          covariate_paths = NULL, output_dir, ...)
    """
    selection = _select_demo_files(base_dir)
    exposures = selection["exposures"]
    mediators = selection["mediators"]
    outcomes = selection["outcomes"]
    covariates = selection["covariates"]
    batches: List[BatchSpec] = []

    if not exposures or not mediators or not outcomes:
        return batches

    demo_out_dir = os.path.join(base_dir, "results_trial", "demo_test")
    for ei, exposure in enumerate(exposures):
        for oi, outcome in enumerate(outcomes):
            mediator = mediators[(ei + oi) % len(mediators)]
            name = f"demo_e{ei+1:02d}_o{oi+1:02d}"
            batches.append(
                BatchSpec(
                    name=name,
                    parameters={
                        "exposure_path": exposure,
                        "mediator_path": mediator,
                        "outcome_path": outcome,
                        "covariate_paths": covariates,
                        # Add any additional keyword args expected by the R function
                    },
                    output_dir=os.path.join(demo_out_dir, name),
                    env={},  # Optional per-batch env vars
                )
            )
    return batches


def _build_full_batches(base_dir: str) -> List[BatchSpec]:
    """
    Build batches for ALL files in the base directory.

    Creates batches for all combinations of:
      - All exposure files from all exposure folders
      - All mediator files
      - All outcome files
      - All covariates

    This is the full processing mode (--full flag).

    Assumes the R entrypoint signature:
        run_main_analysis(exposure_path, mediator_path, outcome_path,
                          covariate_paths = NULL, output_dir, ...)
    """
    expo_dirs = [
        os.path.join(base_dir, "Standardized Circulating human plasma proteome_Data"),
        os.path.join(base_dir, "Standardized Circulating metabolic biomarkers_Data"),
        os.path.join(base_dir, "Circulating inflammatory proteins_Data"),
    ]
    medi_dir = os.path.join(base_dir, "Cerebrospinal fluid metabolomics_Data")
    out_dir = os.path.join(base_dir, "Outcomes")
    covar_dir = os.path.join(base_dir, "Covariates_SES")

    # Collect ALL files from each directory
    exposures: List[str] = []
    for d in expo_dirs:
        exposures.extend(_list_gwas_files(d))

    mediators = _list_gwas_files(medi_dir)
    outcomes = _list_gwas_files(out_dir)
    covariates = _list_gwas_files(covar_dir)

    batches: List[BatchSpec] = []
    if not exposures or not mediators or not outcomes:
        return batches

    # Create output directory for full run
    full_out_dir = os.path.join(base_dir, "results_trial", "full_run")

    # Create batches for all combinations
    batch_counter = 0
    for ei, exposure in enumerate(exposures):
        for mi, mediator in enumerate(mediators):
            for oi, outcome in enumerate(outcomes):
                batch_counter += 1
                name = f"full_e{ei+1:03d}_m{mi+1:03d}_o{oi+1:03d}"
                batches.append(
                    BatchSpec(
                        name=name,
                        parameters={
                            "exposure_path": exposure,
                            "mediator_path": mediator,
                            "outcome_path": outcome,
                            "covariate_paths": covariates,
                        },
                        output_dir=os.path.join(full_out_dir, name),
                        env={},
                    )
                )

    return batches


def _build_optimized_batches(base_dir: str) -> List[BatchSpec]:
    """
    Optimized batch builder with the SAME analytical coverage as full mode
    (all E×M×O combinations) but far fewer batches.

    Strategy:
    - Each batch is an Exposure × Outcome pair.
    - All mediators are passed in a single call, so R loops over every mediator
      internally and still covers the full E×M×O grid.
    - Outputs are grouped under results_trial/optimized_run for separation.
    """
    expo_dirs = [
        os.path.join(base_dir, "Standardized Circulating human plasma proteome_Data"),
        os.path.join(base_dir, "Standardized Circulating metabolic biomarkers_Data"),
        os.path.join(base_dir, "Circulating inflammatory proteins_Data"),
    ]
    medi_dir = os.path.join(base_dir, "Cerebrospinal fluid metabolomics_Data")
    out_dir = os.path.join(base_dir, "Outcomes")
    covar_dir = os.path.join(base_dir, "Covariates_SES")

    # Collect ALL files (no limits)
    exposures: List[str] = []
    for d in expo_dirs:
        exposures.extend(_list_gwas_files(d))

    mediators = _list_gwas_files(medi_dir)
    outcomes = _list_gwas_files(out_dir)
    covariates = _list_gwas_files(covar_dir)

    batches: List[BatchSpec] = []
    if not exposures or not outcomes or not mediators:
        return batches

    # Create output directory
    optimized_out_dir = os.path.join(base_dir, "results_trial", "optimized_run")

    # NOTE: Keep the same E×M×O coverage as full mode but reduce batch count to E×O
    batch_counter = 0
    for oi, outcome in enumerate(outcomes):
        for ei, exposure in enumerate(exposures):
            batch_counter += 1
            name = f"opt_e{ei+1:04d}_o{oi+1:02d}"

            batches.append(
                BatchSpec(
                    name=name,
                    parameters={
                        "exposure_path": exposure,
                        # Pass ALL mediators; R will iterate and produce the full
                        # set of mediator combinations for this exposure/outcome.
                        "mediator_path": mediators,
                        "outcome_path": outcome,
                        "covariate_paths": covariates,
                    },
                    output_dir=os.path.join(optimized_out_dir, name),
                    env={},
                )
            )

    # Metrics (coverage matches full mode; batch count reduced to E×O)
    total_combos = len(exposures) * len(mediators) * len(outcomes)
    total_batches = len(exposures) * len(outcomes)
    print(f"\n[OPTIMIZED MODE] Batch construction complete (full coverage)")
    print(f"  Exposures:          {len(exposures)}")
    print(f"  Mediators (per batch): {len(mediators)}")
    print(f"  Outcomes:           {len(outcomes)}")
    print(f"  Covariates:         {len(covariates)}")
    print(f"  Batches (E×O):      {total_batches:,}")
    print(f"  Combinations covered (E×M×O): {total_combos:,}")
    print(f"  Output dir:         {optimized_out_dir}")
    print("  Note: each batch passes all mediators to R, reducing R session count while retaining full analytical coverage.\n")

    return batches


def _missing_paths(batches: Iterable[BatchSpec]) -> List[str]:
    """
    Collect missing file paths from batch parameters (string-valued).
    """
    missing: List[str] = []
    for batch in batches:
        for value in batch.parameters.values():
            if isinstance(value, str):
                if not os.path.exists(value):
                    missing.append(value)
            elif isinstance(value, Sequence) and not isinstance(value, (bytes, str)):
                for v in value:
                    if isinstance(v, str) and not os.path.exists(v):
                        missing.append(v)
    return missing


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run R batches in parallel using rpy2.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--script",
        default="Main analysis.R",
        help="Path to the R script to source inside each worker.",
    )
    parser.add_argument(
        "--entrypoint",
        default="run_main_analysis",
        help="Name of the R function to call after sourcing the script.",
    )
    parser.add_argument(
        "--processes",
        type=int,
        default=None,
        help="Number of worker processes (defaults to cpu count).",
    )
    parser.add_argument(
        "--example-base-dir",
        type=str,
        default=DEFAULT_EXAMPLE_BASE_DIR,
        help="Run an example batch set constructed from this base dir using the Demo_Test_Analysis.R selection logic. Default mirrors the R config: D:/Projects_data&code/MR_pipeline_demo (outputs to results_trial/demo_test).",
    )
    parser.add_argument(
        "--full",
        action="store_true",
        default=False,
        help="Process all files in the base directory instead of using the demo subset. When enabled, creates batches for all exposures × mediators × outcomes combinations.",
    )
    parser.add_argument(
        "--optimized",
        action="store_true",
        default=False,
        help="Use OPTIMIZED mode with algorithmic improvements. Creates E×O batches instead of E×M×O (passes all mediators per batch) for large reduction in R launches while keeping full coverage.",
    )
    parser.add_argument(
        "--skip-completed",
        action="store_true",
        default=False,
        help="Skip batches that have already been completed (checks batch output markers/files and parent-level completion manifests/logs). Useful for resuming interrupted runs, including after compacting finished batch directories.",
    )
    parser.add_argument(
        "--resume",
        action="store_true",
        default=False,
        help="Alias for --skip-completed. Resume from where the previous run stopped.",
    )
    parser.add_argument(
        "--r-workers",
        type=int,
        default=None,
        help="Override MR_PARALLEL_WORKERS (R-level parallelism per process).",
    )
    parser.add_argument(
        "--dt-threads",
        type=int,
        default=None,
        help="Override MR_DT_THREADS (data.table threads per R process).",
    )
    parser.add_argument(
        "--aggressive",
        action="store_true",
        default=False,
        help="Enable aggressive performance optimizations (larger batches, more caching, higher parallelism).",
    )
    parser.add_argument(
        "--ultra",
        action="store_true",
        default=False,
        help="Enable ULTRA mode: skip non-essential analyses (sensitivity tests, leave-one-out, etc.) for maximum speed.",
    )
    parser.add_argument(
        "--disk-cache",
        action="store_true",
        default=True,  # 默认启用磁盘缓存以提升性能
        help="Enable persistent disk cache for read_gwas (MR_DISK_CACHE=true). Default: enabled for better performance.",
    )
    parser.add_argument(
        "--exposure-cache",
        action="store_true",
        default=True,  # 默认启用exposure缓存以提升性能
        help="Enable persistent exposure cache for preload_exposures (MR_EXPOSURE_CACHE=true). Default: enabled for better performance.",
    )
    parser.add_argument(
        "--exposure-cache-full",
        action="store_true",
        default=False,
        help="Allow caching full exposure GWAS in exposure cache (MR_EXPOSURE_CACHE_FULL=true).",
    )
    parser.add_argument(
        "--no-reverse-ivs-cache",
        action="store_true",
        default=False,
        help="Disable disk cache for reverse MR IVs (MR_REVERSE_IVS_DISK_CACHE=false).",
    )
    parser.add_argument(
        "--skip-reverse-mr",
        action="store_true",
        default=False,
        help="Skip reverse MR bidirectionality test in mediation (MR_SKIP_REVERSE_MR=true).",
    )
    parser.add_argument(
        "--skip-mrlap",
        action="store_true",
        default=False,
        help="Skip MRlap sample-overlap correction in UVMR (MR_RUN_MRLAP=false).",
    )
    parser.add_argument(
        "--monitor-refresh",
        type=float,
        default=0.5,
        help="Progress monitor refresh interval in seconds.",
    )
    parser.add_argument(
        "--monitor-io",
        type=float,
        default=None,
        help="Progress/error log polling interval in seconds (defaults to --monitor-refresh).",
    )
    parser.add_argument(
        "--no-monitor",
        action="store_true",
        default=False,
        help="Disable the live progress monitor to reduce overhead.",
    )
    parser.add_argument(
        "--no-progress-log",
        action="store_true",
        default=False,
        help="Disable writing progress.log from R to reduce I/O overhead.",
    )
    return parser.parse_args()


def main() -> None:
    args = _parse_args()

    # Choose batch building strategy based on flags
    if args.optimized:
        print("Running in OPTIMIZED mode: algorithmic improvements enabled...")
        print("  Batching: Exposure × Outcome (all mediators per batch) to keep full coverage with fewer R sessions.")
        batches = _build_optimized_batches(args.example_base_dir)
        mode_name = "optimized"
    elif args.full:
        print("Running in FULL mode: processing all files in the directory...")
        print("  WARNING: This creates E×M×O batches (may be very large!)")
        batches = _build_full_batches(args.example_base_dir)
        mode_name = "full"
    else:
        print("Running in DEMO mode: processing a subset of files for testing...")
        batches = _build_example_batches(args.example_base_dir)
        mode_name = "demo"

    if not batches:
        print(f"No {mode_name} batches could be constructed. Check that exposure/mediator/outcome files exist under the expected folders.")
        print(f"Base dir tried: {args.example_base_dir}")
        return

    print(f"Total batches to process: {len(batches)}")

    missing = _missing_paths(batches)
    if missing:
        print("Example batch paths not found. Please replace placeholders with real files:")
        for path in missing:
            print(f"  - {path}")
        print("\nProvide your actual base dir with --example-base-dir or edit _build_example_batches().")
        return

    # Determine effective process count (used for auto-tuning R workers/threads)
    skip_completed = args.skip_completed or args.resume
    effective_processes = _effective_processes(args.processes, len(batches))
    cpu_count = os.cpu_count() or 1
    per_process = max(1, cpu_count // effective_processes)

    env_mr_workers = _env_int("MR_PARALLEL_WORKERS")
    env_dt_threads = _env_int("MR_DT_THREADS")

    # 优化：智能并行配置，避免嵌套并行导致的CPU过度订阅
    # 策略：Python多进程 + R内部串行（当Python进程数>1时）
    #       或 Python单进程 + R内部并行（当Python进程数=1时）
    # 激进模式：更激进的并行配置，压榨所有CPU资源
    if args.r_workers is not None and args.r_workers > 0:
        mr_workers = args.r_workers
    elif env_mr_workers is not None and env_mr_workers > 0:
        mr_workers = env_mr_workers
    else:
        # 当Python使用多进程时，R内部使用较少的worker以避免过度订阅
        if effective_processes > 1:
            # 激进模式：允许更多的R worker
            if args.aggressive:
                mr_workers = min(4, per_process)  # 激进：每个进程最多4个worker
            else:
                mr_workers = min(2, per_process)  # 保守：每个进程最多2个worker
        else:
            # 单进程模式：R可以使用所有CPU
            mr_workers = cpu_count

    if args.dt_threads is not None and args.dt_threads > 0:
        dt_threads = args.dt_threads
    elif env_dt_threads is not None and env_dt_threads > 0:
        dt_threads = env_dt_threads
    else:
        # data.table线程数：当R使用并行时保持为1，否则可以使用多线程
        # 激进模式：允许更多的data.table线程
        if args.aggressive:
            dt_threads = 2 if mr_workers > 1 else min(8, per_process)  # 激进：更多线程
        else:
            dt_threads = 1 if mr_workers > 1 else min(4, per_process)  # 保守：标准配置

    mr_workers = max(1, min(cpu_count, int(mr_workers)))
    dt_threads = max(1, min(cpu_count, int(dt_threads)))

    r_parallel_enable = _env_bool("MR_PARALLEL_ENABLE", default=(mr_workers > 1))

    # 激进模式提示
    mode_str = "AGGRESSIVE" if args.aggressive else "STANDARD"
    print(
        f"[TUNING] Mode={mode_str} | processes={effective_processes} | "
        f"R workers={mr_workers} | data.table threads={dt_threads} | "
        f"R parallel={'on' if r_parallel_enable else 'off'}"
    )

    r_env = {
        "MR_SKIP_AUTORUN": "true",
        "BASE_DIR": args.example_base_dir,
        "MR_PARALLEL_WORKERS": str(mr_workers),
        "MR_PARALLEL_ENABLE": "true" if r_parallel_enable else "false",
        "MR_DT_THREADS": str(dt_threads),
    }

    # Optional persistent caches to reduce repeated I/O across batches/processes
    # 优化：默认启用缓存以提升性能，除非用户明确禁用
    if args.disk_cache:
        r_env["MR_DISK_CACHE"] = "true"
        print(f"[CACHE] Disk cache enabled (persistent GWAS read cache)")
    if args.exposure_cache:
        r_env["MR_EXPOSURE_CACHE"] = "true"
        print(f"[CACHE] Exposure cache enabled (persistent exposure preload cache)")
    if args.exposure_cache_full:
        r_env["MR_EXPOSURE_CACHE_FULL"] = "true"
        print(f"[CACHE] Full exposure cache enabled (includes full GWAS data)")
    reverse_cache = not args.no_reverse_ivs_cache
    if reverse_cache:
        r_env["MR_REVERSE_IVS_DISK_CACHE"] = "true"
        print(f"[CACHE] Reverse IVs disk cache enabled (bidirectional MR reuse)")

    if args.skip_reverse_mr:
        r_env["MR_SKIP_REVERSE_MR"] = "true"
        print("[OPT] Reverse MR skipped (MR_SKIP_REVERSE_MR=true)")
    if args.skip_mrlap:
        r_env["MR_RUN_MRLAP"] = "false"
        print("[OPT] MRlap skipped (MR_RUN_MRLAP=false)")

    # 激进模式：设置更大的批量大小
    if args.aggressive:
        r_env["MR_CLUMP_BATCH_SIZE"] = "100"
        print(f"[AGGRESSIVE] Clump batch size set to 100 (default: 50)")

    # ULTRA模式：跳过非关键分析以获得极致速度
    if args.ultra:
        r_env["MR_SKIP_HETEROGENEITY"] = "true"
        r_env["MR_SKIP_PLEIOTROPY"] = "true"
        r_env["MR_SKIP_LOO"] = "true"
        r_env["MR_FAST_MODE"] = "true"
        print(f"[ULTRA] Fast mode enabled - skipping sensitivity tests for maximum speed")
    for var in ("PLINK_BIN", "EUR_BFILE"):
        if os.environ.get(var):
            r_env[var] = os.environ[var]

    config = RRuntimeConfig(
        script_path=args.script,
        entrypoint=args.entrypoint,
        r_options={
            "stringsAsFactors": False,
            # Use multicore futures inside R for lower overhead than multisession
            "mr_future_strategy": "multicore",
            "mr_parallel_enable": r_parallel_enable,
            "mr_parallel_workers": mr_workers,
            "mr_dt_threads": dt_threads,
        },
        env=r_env,
        # Configure R worker count and data.table threads
        dt_threads=dt_threads,
        r_workers=mr_workers,
    )

    monitor_refresh = max(0.1, float(args.monitor_refresh))
    monitor_io = None if args.monitor_io is None else max(0.1, float(args.monitor_io))
    progress_log = not args.no_progress_log
    results = run_batches(
        batches,
        config=config,
        processes=args.processes,
        skip_completed=skip_completed,
        monitor=not args.no_monitor,
        monitor_refresh=monitor_refresh,
        monitor_io=monitor_io,
        progress_log=progress_log,
    )
    failures = [r for r in results if r["status"] == "error"]
    print("\n--- Run summary ---")
    print(f"Total batches: {len(results)}, failures: {len(failures)}")
    for res in results:
        status = res["status"]
        name = res["batch"]
        duration = f"{res.get('duration', 0):.1f}s"
        print(f"{name}: {status} ({duration})")
        if status == "error":
            print(f"  Error: {res['error']}")


if __name__ == "__main__":
    main()
