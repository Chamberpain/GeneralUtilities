#!/usr/bin/env python3
"""
Safely redownload a HYCOM archive with several terminal workers.

Workflow
--------
1. Split the master manifest:
       python parallel_hycom_redownload.py prepare --workers 4

2. Run one worker in each terminal:
       python parallel_hycom_redownload.py worker --workers 4 --index 0
       python parallel_hycom_redownload.py worker --workers 4 --index 1
       python parallel_hycom_redownload.py worker --workers 4 --index 2
       python parallel_hycom_redownload.py worker --workers 4 --index 3

3. Merge the worker manifests:
       python parallel_hycom_redownload.py merge --workers 4

Legacy CEC entrypoints remain supported. Pass --module-dir to select an
alternative hycom_download.py repair implementation.
"""

from __future__ import annotations

import argparse
import contextlib
import os
import shutil
import sys
import importlib.util
import tempfile
from pathlib import Path
from typing import Iterator

import pandas as pd


DEFAULT_HYCOM_DIR = Path(
    "/Users/paulchamberlain/Projects/GeneralUtilities/"
    "DataDir/CEC/Data/hycom_currents_humboldt"
)
DEFAULT_DATA_ROOT = Path(
    "/Users/paulchamberlain/Projects/GeneralUtilities/DataDir/CEC/Data"
)

ROW_ID_COLUMN = "_parallel_master_row"


def atomic_csv_write(frame: pd.DataFrame, destination: Path) -> None:
    """Write a CSV atomically in the destination directory."""
    destination.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary_name = tempfile.mkstemp(
        prefix=destination.name + ".",
        suffix=".tmp",
        dir=destination.parent,
    )
    try:
        with os.fdopen(fd, "w", encoding="utf-8", newline="") as handle:
            frame.to_csv(handle, index=False)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary_name, destination)
    finally:
        Path(temporary_name).unlink(missing_ok=True)


def worker_manifest_path(hycom_dir: Path, index: int) -> Path:
    return hycom_dir / f"download_manifest.worker_{index}.csv"


def resolve_output_path(raw_path: str, data_root: Path) -> Path:
    path = Path(raw_path).expanduser()
    if path.is_absolute():
        return path.resolve()
    return (data_root / path).resolve()


def ensure_beneath(path: Path, directory: Path) -> None:
    try:
        path.relative_to(directory)
    except ValueError as exc:
        raise ValueError(
            f"Refusing to touch a manifest path outside {directory}: {path}"
        ) from exc


@contextlib.contextmanager
def exclusive_target_lock(target: Path) -> Iterator[None]:
    """
    Lock one target for the duration of its replacement.

    This is an additional guard against accidentally starting the same worker
    twice. It does not serialize different files.
    """
    try:
        import fcntl
    except ImportError as exc:  # pragma: no cover - intended for macOS/Linux
        raise RuntimeError(
            "This helper requires the Unix fcntl module (available on macOS/Linux)."
        ) from exc

    lock_path = target.with_name(target.name + ".parallel.lock")
    lock_path.parent.mkdir(parents=True, exist_ok=True)

    with open(lock_path, "a+", encoding="utf-8") as lock_handle:
        fcntl.flock(lock_handle.fileno(), fcntl.LOCK_EX)
        try:
            yield
        finally:
            fcntl.flock(lock_handle.fileno(), fcntl.LOCK_UN)


def prepare(args: argparse.Namespace) -> int:
    hycom_dir = args.hycom_dir.expanduser().resolve()
    manifest = args.manifest.expanduser().resolve()
    workers = args.workers

    if workers < 1:
        raise ValueError("--workers must be at least 1")
    if not manifest.is_file():
        raise FileNotFoundError(f"Master manifest not found: {manifest}")

    frame = pd.read_csv(manifest, dtype=str, keep_default_na=False)
    required = {"date", "output_path"}
    missing = required - set(frame.columns)
    if missing:
        raise ValueError(
            f"Master manifest is missing required columns: {sorted(missing)}"
        )
    if ROW_ID_COLUMN in frame.columns and not args.force:
        raise ValueError(
            f"{ROW_ID_COLUMN!r} already exists. Use --force only after checking "
            "that no prior parallel download is active."
        )

    backup = manifest.with_name(manifest.stem + ".before_parallel.csv")
    destinations = [worker_manifest_path(hycom_dir, index) for index in range(workers)]
    # Check the entire run before creating a backup or any worker manifests.
    # Otherwise a collision in a later worker leaves an unusable partial split.
    if backup.exists() and not args.force:
        raise FileExistsError(
            f"Backup already exists: {backup}\n"
            "This usually means a prior parallel run has not been merged or cleaned up."
        )
    for destination in destinations:
        if destination.exists() and not args.force:
            raise FileExistsError(
                f"Worker manifest already exists: {destination}. "
                "Use --force only when no worker is running."
            )
        if destination.is_dir():
            raise IsADirectoryError(
                f"Worker manifest path is a directory: {destination}"
            )
    if backup.is_dir():
        raise IsADirectoryError(f"Backup path is a directory: {backup}")
    if manifest in destinations:
        raise ValueError("Master manifest must not also be a worker manifest.")

    # Remove any stale helper column before assigning a fresh stable row number.
    frame = frame.drop(columns=[ROW_ID_COLUMN], errors="ignore")
    frame.insert(0, ROW_ID_COLUMN, [str(i) for i in range(len(frame))])

    # Worker manifests represent a new full-redownload run. Reset only the
    # completion markers used by this helper; request URLs, dates, checksums,
    # and all other source metadata remain available to repair_missing_files.
    frame["status"] = "parallel_pending"
    frame["repair_status"] = "pending"
    frame["repair_error"] = ""

    shutil.copy2(manifest, backup)

    for index, destination in enumerate(destinations):
        part = frame.iloc[index::workers]
        atomic_csv_write(part, destination)
        print(f"worker {index}: {len(part)} rows -> {destination}")

    print(f"\nMaster manifest left unchanged: {manifest}")
    print(f"Backup created: {backup}")
    return 0


def _load_repair_module(module_dir: Path):
    """Load a selected repair module without consulting Python's import cache."""
    module_dir = module_dir.expanduser().resolve()
    if module_dir == Path(__file__).resolve().parent:
        from . import cec_hycom_download
        return cec_hycom_download
    source_path = module_dir / "hycom_download.py"
    if not source_path.is_file():
        raise FileNotFoundError(
            f"hycom_download.py was not found in --module-dir {module_dir}"
        )
    spec = importlib.util.spec_from_file_location("hycom_download", source_path)
    module = importlib.util.module_from_spec(spec)
    missing = object()
    previous = sys.modules.get(spec.name, missing)
    try:
        spec.loader.exec_module(module)
        # CEC's compatibility entrypoint aliases the canonical shared module.
        return getattr(module, "_implementation", module)
    finally:
        if previous is missing:
            sys.modules.pop(spec.name, None)
        else:
            sys.modules[spec.name] = previous


def import_repair_function(module_dir: Path):
    """Return the repair function from the explicitly selected module directory."""
    return _load_repair_module(module_dir).repair_missing_files


def _run_worker(args: argparse.Namespace, repair_factory, *, on_failure=None) -> int:
    """Process a worker manifest with shared locking, backup and recovery rules.

    repair_factory is resolved after manifest/index preflight and returns a
    callable accepting target, manifest, data root and original-file template.
    The callable may return a strategy label for the completion message.
    """
    hycom_dir = args.hycom_dir.expanduser().resolve()
    data_root = args.data_root.expanduser().resolve()
    manifest = worker_manifest_path(hycom_dir, args.index)

    if args.index < 0 or args.index >= args.workers:
        raise ValueError("--index must be between 0 and --workers - 1")
    if not manifest.is_file():
        raise FileNotFoundError(
            f"Worker manifest not found: {manifest}\nRun the prepare subcommand first."
        )

    repair_target = repair_factory()

    frame = pd.read_csv(manifest, dtype=str, keep_default_na=False)
    required = {"date", "output_path", ROW_ID_COLUMN}
    missing = required - set(frame.columns)
    if missing:
        raise ValueError(
            f"Worker manifest is missing required columns: {sorted(missing)}"
        )

    failures: list[tuple[str, str]] = []
    skipped = 0
    completed = 0

    for sequence, (_, row) in enumerate(frame.iterrows(), start=1):
        target = resolve_output_path(row["output_path"], data_root)
        ensure_beneath(target, hycom_dir)

        print(
            f"\n[worker {args.index}: {sequence}/{len(frame)}] {target.name}",
            flush=True,
        )

        with exclusive_target_lock(target):
            # A restarted worker skips rows that it previously completed and
            # whose replacement file still exists.
            status = str(row.get("status", "")).strip().lower()
            repair_status = str(row.get("repair_status", "")).strip().lower()
            if (
                not args.force_completed
                and status == "downloaded"
                and repair_status == "repaired"
                and target.is_file()
                and target.stat().st_size > 0
            ):
                print("Already completed by this worker; skipping.", flush=True)
                skipped += 1
                continue

            target.parent.mkdir(parents=True, exist_ok=True)
            backup = target.with_name(
                target.stem + ".before-redownload" + target.suffix
            )

            # Recover cleanly from an interrupted prior attempt.
            if backup.is_file() and not target.exists():
                template = backup
            elif backup.exists() and target.exists():
                raise FileExistsError(
                    f"Both target and backup exist:\n  {target}\n  {backup}\n"
                    "Inspect them before restarting this worker."
                )
            elif target.is_file() and target.stat().st_size > 0:
                os.replace(target, backup)
                template = backup
            elif target.exists():
                if target.is_file():
                    target.unlink()
                    template = backup if backup.is_file() else None
                else:
                    raise ValueError(
                        f"Manifest destination is not a regular file: {target}"
                    )
            else:
                template = backup if backup.is_file() else None

            try:
                mode = repair_target(target, manifest, data_root, template)
            except Exception as exc:
                # Restore the original data file. The worker manifest retains
                # the failed status written by repair_missing_files.
                if backup.is_file():
                    if target.exists():
                        target.unlink()
                    os.replace(backup, target)

                failures.append((str(target), f"{type(exc).__name__}: {exc}"))
                print(f"FAILED: {target}\n{exc}", flush=True)

                if on_failure is not None:
                    on_failure()
                if not args.continue_on_error:
                    break
            else:
                if backup.exists():
                    backup.unlink()
                completed += 1
                label = f" ({mode})" if mode else ""
                print(f"COMPLETE{label}: {target}", flush=True)

    print(
        f"\nworker {args.index} summary: "
        f"completed={completed}, skipped={skipped}, failed={len(failures)}"
    )
    if failures:
        for path, message in failures:
            print(f"  {path}: {message}")
        return 1
    return 0


def worker(args: argparse.Namespace) -> int:
    def repair_factory():
        repair_missing_files = import_repair_function(args.module_dir)

        def repair_target(target, manifest, data_root, template):
            repair_missing_files(
                [target],
                manifest_path=manifest,
                data_root=data_root,
                template_path=template,
                retries=args.retries,
                timeout=(args.connect_timeout, args.read_timeout),
                backoff_seconds=args.backoff,
            )

        return repair_target

    return _run_worker(args, repair_factory)


def merge(args: argparse.Namespace) -> int:
    hycom_dir = args.hycom_dir.expanduser().resolve()
    data_root = args.data_root.expanduser().resolve()
    master = args.manifest.expanduser().resolve()

    if args.workers < 1:
        raise ValueError("--workers must be at least 1")
    if not master.is_file():
        raise FileNotFoundError(f"Master manifest not found: {master}")

    worker_paths = [
        worker_manifest_path(hycom_dir, index) for index in range(args.workers)
    ]
    if master in worker_paths:
        raise ValueError("Master manifest must not also be a worker manifest.")
    missing_paths = [path for path in worker_paths if not path.is_file()]
    if missing_paths:
        raise FileNotFoundError(
            "Missing worker manifests:\n" + "\n".join(map(str, missing_paths))
        )

    parts = [
        pd.read_csv(path, dtype=str, keep_default_na=False) for path in worker_paths
    ]
    merged = pd.concat(parts, ignore_index=True, sort=False)

    original = pd.read_csv(master, dtype=str, keep_default_na=False)
    required = {"date", "output_path"}
    if not required.issubset(original.columns):
        raise ValueError("Master manifest requires date and output_path columns")
    missing = required.union({ROW_ID_COLUMN}) - set(merged.columns)
    if missing:
        raise ValueError(
            f"Worker manifests are missing required columns: {sorted(missing)}"
        )

    # These IDs refer to positions in the original master, not arbitrary sort
    # keys. Require every position exactly once, including for an empty run.
    row_ids = pd.to_numeric(merged[ROW_ID_COLUMN], errors="raise")
    if row_ids.sort_values().tolist() != list(range(len(original))):
        raise ValueError(
            "Worker manifests must contain each master row identifier exactly once."
        )
    merged = merged.iloc[row_ids.argsort()].reset_index(drop=True)

    original_keys = list(zip(original["date"], original["output_path"]))
    merged_keys = list(zip(merged["date"], merged["output_path"]))
    if original_keys != merged_keys:
        raise ValueError(
            "Worker manifests do not contain exactly the same date/output_path "
            "records as the master manifest."
        )

    status = merged.get("status", pd.Series("", index=merged.index)).str.lower()
    repair_status = merged.get(
        "repair_status", pd.Series("", index=merged.index)
    ).str.lower()
    complete_mask = status.eq("downloaded") & repair_status.eq("repaired")
    failure_count = int(status.eq("download_failed").sum())
    incomplete_count = int((~complete_mask).sum())

    missing_outputs = []
    for raw_path in merged["output_path"]:
        output = resolve_output_path(raw_path, data_root)
        ensure_beneath(output, hycom_dir)
        if not output.is_file() or output.stat().st_size == 0:
            missing_outputs.append(str(output))

    if (incomplete_count or missing_outputs) and not args.allow_incomplete:
        details = [
            f"{incomplete_count} worker rows are not marked downloaded/repaired",
            f"{failure_count} rows are marked download_failed",
            f"{len(missing_outputs)} output files are missing or empty",
        ]
        raise RuntimeError(
            "; ".join(details)
            + ". Rerun the affected workers, or pass --allow-incomplete "
            "to merge the partial state intentionally."
        )

    merged = merged.drop(columns=[ROW_ID_COLUMN])
    premerge_backup = master.with_name(master.stem + ".before_merge.csv")
    if not premerge_backup.exists():
        shutil.copy2(master, premerge_backup)

    atomic_csv_write(merged, master)
    print(f"Merged {len(merged)} rows into {master}")
    print(f"Rows not marked downloaded/repaired: {incomplete_count}")
    print(f"Rows marked download_failed: {failure_count}")
    print(f"Missing or empty output files: {len(missing_outputs)}")
    print(f"Pre-merge backup: {premerge_backup}")

    if args.cleanup:
        for path in worker_paths:
            path.unlink(missing_ok=True)
            path.with_suffix(path.suffix + ".bak").unlink(missing_ok=True)
        print("Worker manifests and their .bak files were removed.")

    return 0


def _build_parser(
    *,
    worker_func,
    description,
    hycom_dir,
    data_root,
    retries=5,
    connect_timeout=15,
    read_timeout=180,
    backoff=2,
    worker_help=None,
):
    """Create common archive commands and return the worker parser to extend."""
    parser = argparse.ArgumentParser(description=description)
    parser.set_defaults(func=None)

    common = argparse.ArgumentParser(add_help=False)
    common.add_argument(
        "--hycom-dir",
        type=Path,
        default=hycom_dir,
        help=f"HYCOM output directory (default: {hycom_dir})",
    )
    common.add_argument(
        "--data-root",
        type=Path,
        default=data_root,
        help=f"Base for relative output_path values (default: {data_root})",
    )
    common.add_argument(
        "--manifest",
        type=Path,
        default=hycom_dir / "download_manifest.csv",
        help="Master download_manifest.csv",
    )
    common.add_argument(
        "--workers",
        type=int,
        default=4,
        help="Number of terminal workers (default: 4)",
    )

    subparsers = parser.add_subparsers(dest="command")

    prepare_parser = subparsers.add_parser(
        "prepare",
        parents=[common],
        help="Create one disjoint worker manifest per terminal.",
    )
    prepare_parser.add_argument(
        "--force",
        action="store_true",
        help="Replace existing worker manifests and the parallel backup.",
    )
    prepare_parser.set_defaults(func=prepare)

    worker_parser = subparsers.add_parser(
        "worker",
        parents=[common],
        help="Run one terminal worker.",
    )
    worker_parser.add_argument("--index", type=int, required=True)
    worker_parser.add_argument(
        "--module-dir",
        type=Path,
        default=Path.cwd(),
        help="Directory containing hycom_download.py (default: current directory)",
    )
    worker_help = worker_help or {}
    for flag, value_type, default in (
        ("--retries", int, retries),
        ("--connect-timeout", float, connect_timeout),
        ("--read-timeout", float, read_timeout),
        ("--backoff", float, backoff),
    ):
        worker_parser.add_argument(
            flag, type=value_type, default=default, help=worker_help.get(flag)
        )
    worker_parser.add_argument(
        "--continue-on-error",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Continue to later dates after one date fails (default: true)",
    )
    worker_parser.add_argument(
        "--force-completed",
        action="store_true",
        help="Redownload rows already marked downloaded/repaired by this worker.",
    )
    worker_parser.set_defaults(func=worker_func)

    merge_parser = subparsers.add_parser(
        "merge",
        parents=[common],
        help="Merge completed worker manifests into the master manifest.",
    )
    merge_parser.add_argument(
        "--allow-incomplete",
        action="store_true",
        help=(
            "Merge even if rows are unfinished or outputs are missing. "
            "Normally the merge refuses a partial run."
        ),
    )
    merge_parser.add_argument(
        "--cleanup",
        action="store_true",
        help="Remove worker manifests after a successful merge.",
    )
    merge_parser.set_defaults(func=merge)

    return parser, worker_parser


def build_parser() -> argparse.ArgumentParser:
    parser, _ = _build_parser(
        worker_func=worker,
        description=(
            "Split a HYCOM manifest among independent terminal workers, "
            "redownload each worker's files, and merge the results safely."
        ),
        hycom_dir=DEFAULT_HYCOM_DIR,
        data_root=DEFAULT_DATA_ROOT,
    )
    return parser


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()
    if args.func is None:
        parser.print_help()
        return 2
    return int(args.func(args))


if __name__ == "__main__":
    raise SystemExit(main())
