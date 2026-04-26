#!/usr/bin/env python3
"""
Compact optimized_run batch outputs by:
1. Merging UVMR / MVMR / mediation result CSVs into full_* CSVs.
2. Recording all completed batch names into completed_batches.txt.
3. Optionally deleting completed batch directories afterward.

Completion detection matches parallel_r_runner.py so that
--skip-completed / --resume still work after compaction.
"""

from __future__ import annotations

import argparse
import shutil
import sys
from collections import OrderedDict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Optional

try:
    import pandas as pd
except ImportError as exc:  # pragma: no cover - runtime dependency guard
    raise SystemExit(
        "pandas is required to run this script. Install it first, then retry."
    ) from exc

from parallel_r_runner import _is_batch_completed


@dataclass(frozen=True)
class ResultSpec:
    key: str
    log_name: str
    result_name: str
    merged_name: str


@dataclass(frozen=True)
class LogEntry:
    spec_key: str
    line_number: int
    csv_path: Path
    batch_dir: Path
    batch_name: str
    rows: Optional[int]


RESULT_SPECS = (
    ResultSpec("uvmr", "uvmr.log", "uvmr_results.csv", "full_uvmr_results.csv"),
    ResultSpec("mvmr", "mvmr.log", "mvmr_results.csv", "full_mvmr_results.csv"),
    ResultSpec("mediation", "mediation.log", "mediation_results.csv", "full_mediation_results.csv"),
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Merge optimized_run result CSVs (UVMR/MVMR/mediation), write a "
            "completed batch manifest, and optionally delete completed batch dirs."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "batch_root",
        type=Path,
        help="Root directory containing opt_* batch folders, e.g. optimized_run",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Directory for merged full_* CSVs. Defaults to batch_root.",
    )
    parser.add_argument(
        "--manifest",
        type=Path,
        default=None,
        help="Completed batch manifest path. Defaults to <batch_root>/completed_batches.txt",
    )
    parser.add_argument(
        "--uvmr-log",
        type=Path,
        default=None,
        help="Override path to uvmr.log. Defaults to <batch_root>/uvmr.log",
    )
    parser.add_argument(
        "--mvmr-log",
        type=Path,
        default=None,
        help="Override path to mvmr.log. Defaults to <batch_root>/mvmr.log",
    )
    parser.add_argument(
        "--mediation-log",
        type=Path,
        default=None,
        help="Override path to mediation.log. Defaults to <batch_root>/mediation.log",
    )
    parser.add_argument(
        "--delete-batch-dirs",
        action="store_true",
        help="Delete completed opt_* batch directories after outputs/manifests are written.",
    )
    parser.add_argument(
        "--strict-missing",
        action="store_true",
        help="Fail if a logged CSV is missing instead of skipping it with a warning.",
    )
    parser.add_argument(
        "--add-source-columns",
        action="store_true",
        help="Add source_batch and source_csv columns to merged outputs.",
    )
    return parser.parse_args()


def resolve_log_path(args: argparse.Namespace, spec: ResultSpec, batch_root: Path) -> Optional[Path]:
    override = getattr(args, f"{spec.key}_log")
    if override is not None:
        return override.expanduser().resolve()

    candidate = batch_root / spec.log_name
    if candidate.exists():
        return candidate

    return None


def parse_log_entries(log_path: Path, spec: ResultSpec) -> list[LogEntry]:
    entries_by_csv: "OrderedDict[str, LogEntry]" = OrderedDict()

    with log_path.open("r", encoding="utf-8") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.strip()
            if not line or "\t" not in line:
                continue

            _, message = line.split("\t", 1)
            parts = message.split()
            if len(parts) < 2 or parts[0] != spec.result_name:
                continue

            csv_path = Path(parts[1])
            batch_dir = csv_path.parent
            batch_name = batch_dir.name
            if not batch_name:
                continue

            rows = None
            if len(parts) >= 4 and parts[-2] == "rows":
                try:
                    rows = int(parts[-1])
                except ValueError:
                    rows = None

            entries_by_csv[str(csv_path)] = LogEntry(
                spec_key=spec.key,
                line_number=line_number,
                csv_path=csv_path,
                batch_dir=batch_dir,
                batch_name=batch_name,
                rows=rows,
            )

    return list(entries_by_csv.values())


def list_batch_dirs(batch_root: Path) -> list[Path]:
    if not batch_root.exists() or not batch_root.is_dir():
        return []

    return sorted(
        [
            child
            for child in batch_root.iterdir()
            if child.is_dir() and child.name.startswith("opt_")
        ],
        key=lambda path: path.name,
    )


def collect_completed_batches(
    batch_root: Path,
    entries_by_spec: dict[str, list[LogEntry]],
) -> "OrderedDict[str, Path]":
    completed: "OrderedDict[str, Path]" = OrderedDict()

    for batch_dir in list_batch_dirs(batch_root):
        if _is_batch_completed(str(batch_dir)):
            completed[batch_dir.name] = batch_dir

    for entries in entries_by_spec.values():
        for entry in entries:
            completed.setdefault(entry.batch_name, batch_root / entry.batch_name)

    return completed


def read_frames(
    entries: Iterable[LogEntry],
    *,
    strict_missing: bool,
    add_source_columns: bool,
    batch_root: Path,
) -> tuple[list["pd.DataFrame"], list[LogEntry]]:
    frames: list["pd.DataFrame"] = []
    missing_entries: list[LogEntry] = []

    for entry in entries:
        csv_path = entry.csv_path
        if not csv_path.exists():
            remapped_path = batch_root / entry.batch_name / entry.csv_path.name
            if remapped_path.exists():
                csv_path = remapped_path

        if not csv_path.exists():
            missing_entries.append(entry)
            if strict_missing:
                raise FileNotFoundError(
                    f"Logged CSV is missing: {entry.csv_path} (log line {entry.line_number})"
                )
            continue

        frame = pd.read_csv(csv_path, low_memory=False)
        if add_source_columns:
            frame.insert(0, "source_csv", str(csv_path))
            frame.insert(0, "source_batch", entry.batch_name)
        frames.append(frame)

    return frames, missing_entries


def write_manifest(manifest_path: Path, batch_names: Iterable[str]) -> None:
    manifest_path.parent.mkdir(parents=True, exist_ok=True)

    seen = set()
    ordered_names = []
    for batch_name in batch_names:
        if batch_name in seen:
            continue
        seen.add(batch_name)
        ordered_names.append(batch_name)

    with manifest_path.open("w", encoding="utf-8") as handle:
        handle.write("# One completed batch name per line.\n")
        for batch_name in ordered_names:
            handle.write(f"{batch_name}\n")


def delete_batch_dirs(batch_dirs: Iterable[Path], batch_root: Path) -> int:
    removed = 0
    expected_root = batch_root.resolve()
    seen = set()

    for batch_dir in batch_dirs:
        normalized = str(batch_dir)
        if normalized in seen:
            continue
        seen.add(normalized)

        resolved_batch_dir = batch_dir.resolve()
        if resolved_batch_dir.parent != expected_root:
            print(
                f"[WARN] Refusing to delete unexpected directory outside batch root: {batch_dir}",
                file=sys.stderr,
            )
            continue

        if not resolved_batch_dir.exists():
            continue

        shutil.rmtree(resolved_batch_dir)
        removed += 1

    return removed


def merge_one_result_type(
    spec: ResultSpec,
    entries: list[LogEntry],
    *,
    batch_root: Path,
    output_dir: Path,
    strict_missing: bool,
    add_source_columns: bool,
) -> tuple[int, int, int, Optional[Path]]:
    frames, missing_entries = read_frames(
        entries,
        strict_missing=strict_missing,
        add_source_columns=add_source_columns,
        batch_root=batch_root,
    )

    if missing_entries:
        print(
            f"[WARN] {spec.key}: missing CSV files skipped: {len(missing_entries)}",
            file=sys.stderr,
        )
        for entry in missing_entries[:10]:
            print(f"  - {entry.csv_path}", file=sys.stderr)
        if len(missing_entries) > 10:
            print(f"  ... and {len(missing_entries) - 10} more", file=sys.stderr)

    if not frames:
        return len(entries), 0, 0, None

    merged = pd.concat(frames, ignore_index=True, sort=False)
    output_dir.mkdir(parents=True, exist_ok=True)
    output_path = output_dir / spec.merged_name
    merged.to_csv(output_path, index=False)
    return len(entries), len(frames), len(merged), output_path


def main() -> int:
    args = parse_args()
    batch_root = args.batch_root.expanduser().resolve()
    output_dir = (
        args.output_dir.expanduser().resolve()
        if args.output_dir is not None
        else batch_root
    )
    manifest_path = (
        args.manifest.expanduser().resolve()
        if args.manifest is not None
        else batch_root / "completed_batches.txt"
    )

    entries_by_spec: dict[str, list[LogEntry]] = {}
    for spec in RESULT_SPECS:
        log_path = resolve_log_path(args, spec, batch_root)
        if log_path is None or not log_path.exists():
            entries_by_spec[spec.key] = []
            print(f"[WARN] {spec.key}: log not found, skipping log merge.")
            continue

        entries = parse_log_entries(log_path, spec)
        entries_by_spec[spec.key] = entries
        print(f"[OK] {spec.key}: parsed log entries: {len(entries)}")

    completed_batches = collect_completed_batches(batch_root, entries_by_spec)
    write_manifest(manifest_path, completed_batches.keys())

    for spec in RESULT_SPECS:
        entries = entries_by_spec[spec.key]
        if not entries:
            continue

        log_count, merged_file_count, row_count, output_path = merge_one_result_type(
            spec,
            entries,
            batch_root=batch_root,
            output_dir=output_dir,
            strict_missing=args.strict_missing,
            add_source_columns=args.add_source_columns,
        )
        print(f"[OK] {spec.key}: log entries: {log_count}")
        print(f"[OK] {spec.key}: merged CSV files: {merged_file_count}")
        if output_path is not None:
            print(f"[OK] {spec.key}: rows written: {row_count}")
            print(f"[OK] {spec.key}: output: {output_path}")
        else:
            print(f"[WARN] {spec.key}: no readable CSV files were available; skipped merged CSV write.")

    removed = 0
    if args.delete_batch_dirs:
        removed = delete_batch_dirs(completed_batches.values(), batch_root)

    print(f"[OK] Batch root: {batch_root}")
    print(f"[OK] Completed batches recorded: {len(completed_batches)}")
    print(f"[OK] Manifest: {manifest_path}")
    if args.delete_batch_dirs:
        print(f"[OK] Deleted batch directories: {removed}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
