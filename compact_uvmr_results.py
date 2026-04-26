#!/usr/bin/env python3
"""
Merge per-batch UVMR result CSVs listed in uvmr.log into one full CSV, then
optionally delete the original batch directories after writing lightweight
completion markers for future --skip-completed / --resume runs.
"""

from __future__ import annotations

import argparse
import shutil
import sys
from collections import OrderedDict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, List, Optional

try:
    import pandas as pd
except ImportError as exc:  # pragma: no cover - runtime dependency guard
    raise SystemExit(
        "pandas is required to run this script. Install it first, then retry."
    ) from exc

from parallel_r_runner import _is_batch_completed


@dataclass(frozen=True)
class LogEntry:
    line_number: int
    csv_path: Path
    batch_dir: Path
    batch_name: str
    rows: Optional[int]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Read uvmr.log, merge all logged uvmr_results.csv files into one "
            "full_uvmr_results.csv, and optionally delete the original batch dirs."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "log_path",
        type=Path,
        help="Path to uvmr.log",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="Output path for the merged CSV. Defaults to <log_dir>/full_uvmr_results.csv",
    )
    parser.add_argument(
        "--manifest",
        type=Path,
        default=None,
        help="Completed batch manifest. Defaults to <log_dir>/completed_batches.txt",
    )
    parser.add_argument(
        "--delete-batch-dirs",
        action="store_true",
        help="Delete original batch directories after the merged CSV and manifest are written.",
    )
    parser.add_argument(
        "--strict-missing",
        action="store_true",
        help="Fail if a logged CSV is missing instead of skipping it with a warning.",
    )
    parser.add_argument(
        "--add-source-columns",
        action="store_true",
        help="Add source_batch and source_csv columns to the merged output.",
    )
    parser.add_argument(
        "--batch-root",
        type=Path,
        default=None,
        help=(
            "Root directory containing opt_* batch folders. If provided, the script "
            "scans all completed batch dirs there, including batches missing "
            "uvmr_results.csv but finished according to progress.log/.completed."
        ),
    )
    return parser.parse_args()


def parse_log_entries(log_path: Path) -> List[LogEntry]:
    entries_by_csv: "OrderedDict[str, LogEntry]" = OrderedDict()

    with log_path.open("r", encoding="utf-8") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.strip()
            if not line:
                continue

            if "\t" not in line:
                continue

            _, message = line.split("\t", 1)
            parts = message.split()
            if len(parts) < 2:
                continue
            if parts[0] != "uvmr_results.csv":
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
                line_number=line_number,
                csv_path=csv_path,
                batch_dir=batch_dir,
                batch_name=batch_name,
                rows=rows,
            )

    return list(entries_by_csv.values())


def infer_batch_parent(entries: Iterable[LogEntry]) -> Path:
    """
    Infer the shared parent directory that contains all batch folders.
    """
    parents = {entry.batch_dir.parent for entry in entries}
    if not parents:
        raise ValueError("No batch directories found in log entries.")
    if len(parents) != 1:
        parent_list = ", ".join(str(parent) for parent in sorted(parents))
        raise ValueError(
            "Logged batch directories do not share a single parent directory: "
            f"{parent_list}"
        )
    return next(iter(parents))


def list_batch_dirs(batch_parent: Path) -> List[Path]:
    if not batch_parent.exists() or not batch_parent.is_dir():
        return []

    return sorted(
        [
            child
            for child in batch_parent.iterdir()
            if child.is_dir() and child.name.startswith("opt_")
        ],
        key=lambda path: path.name,
    )


def collect_completed_batches(
    batch_parent: Path,
    entries: Iterable[LogEntry],
) -> "OrderedDict[str, Path]":
    completed: "OrderedDict[str, Path]" = OrderedDict()

    for batch_dir in list_batch_dirs(batch_parent):
        if _is_batch_completed(str(batch_dir)):
            completed[batch_dir.name] = batch_dir

    for entry in entries:
        completed.setdefault(entry.batch_name, entry.batch_dir)

    return completed


def read_frames(
    entries: Iterable[LogEntry],
    *,
    strict_missing: bool,
    add_source_columns: bool,
    batch_root: Optional[Path] = None,
) -> tuple[list["pd.DataFrame"], list[LogEntry]]:
    frames: list["pd.DataFrame"] = []
    missing_entries: list[LogEntry] = []

    for entry in entries:
        csv_path = entry.csv_path
        if not csv_path.exists() and batch_root is not None:
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


def delete_batch_dirs(batch_dirs: Iterable[Path], expected_parent: Path) -> int:
    removed = 0
    expected_parent = expected_parent.resolve()
    seen = set()

    for batch_dir in batch_dirs:
        normalized = str(batch_dir)
        if normalized in seen:
            continue
        seen.add(normalized)

        try:
            resolved_batch_dir = batch_dir.resolve()
        except FileNotFoundError:
            continue

        if resolved_batch_dir.parent != expected_parent:
            print(
                f"[WARN] Refusing to delete unexpected directory outside batch parent: {batch_dir}",
                file=sys.stderr,
            )
            continue

        if not resolved_batch_dir.exists():
            continue

        shutil.rmtree(resolved_batch_dir)
        removed += 1

    return removed


def main() -> int:
    args = parse_args()
    log_path = args.log_path.expanduser().resolve()
    if not log_path.exists():
        raise SystemExit(f"uvmr.log not found: {log_path}")

    log_dir = log_path.parent
    entries = parse_log_entries(log_path)
    if not entries:
        raise SystemExit(f"No uvmr_results.csv entries found in: {log_path}")
    batch_parent = (
        args.batch_root.expanduser().resolve()
        if args.batch_root is not None
        else infer_batch_parent(entries)
    )
    completed_batches = collect_completed_batches(batch_parent, entries)

    output_path = (
        args.output.expanduser().resolve()
        if args.output is not None
        else log_dir / "full_uvmr_results.csv"
    )
    manifest_path = (
        args.manifest.expanduser().resolve()
        if args.manifest is not None
        else batch_parent / "completed_batches.txt"
    )

    frames, missing_entries = read_frames(
        entries,
        strict_missing=args.strict_missing,
        add_source_columns=args.add_source_columns,
        batch_root=batch_parent,
    )
    write_manifest(manifest_path, completed_batches.keys())

    merged = None
    if frames:
        merged = pd.concat(frames, ignore_index=True, sort=False)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        merged.to_csv(output_path, index=False)
    elif not args.delete_batch_dirs:
        raise SystemExit("No readable UVMR CSV files were found; nothing was merged.")

    removed = 0
    if args.delete_batch_dirs:
        removed = delete_batch_dirs(completed_batches.values(), batch_parent)

    print(f"[OK] Parsed log entries: {len(entries)}")
    print(f"[OK] Completed batches recorded: {len(completed_batches)}")
    print(f"[OK] Merged CSV files: {len(frames)}")
    if merged is not None:
        print(f"[OK] Rows written: {len(merged)}")
        print(f"[OK] Output CSV: {output_path}")
    else:
        print("[WARN] No readable UVMR CSV files were available; skipped merged CSV write.")
    print(f"[OK] Batch parent: {batch_parent}")
    print(f"[OK] Manifest: {manifest_path}")
    if missing_entries:
        print(f"[WARN] Missing CSV files skipped: {len(missing_entries)}", file=sys.stderr)
        for entry in missing_entries[:10]:
            print(f"  - {entry.csv_path}", file=sys.stderr)
        if len(missing_entries) > 10:
            print(f"  ... and {len(missing_entries) - 10} more", file=sys.stderr)
    if args.delete_batch_dirs:
        print(f"[OK] Deleted batch directories: {removed}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
