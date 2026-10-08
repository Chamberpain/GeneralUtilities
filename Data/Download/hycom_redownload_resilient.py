#!/usr/bin/env python3
"""
Safely redownload a HYCOM archive with several terminal workers.

Workflow
--------
1. Split the master manifest:
       python parallel_hycom_redownload_resilient.py prepare --workers 4

2. Run one worker in each terminal:
       python parallel_hycom_redownload_resilient.py worker --workers 4 --index 0
       python parallel_hycom_redownload_resilient.py worker --workers 4 --index 1
       python parallel_hycom_redownload_resilient.py worker --workers 4 --index 2
       python parallel_hycom_redownload_resilient.py worker --workers 4 --index 3

3. Merge the worker manifests:
       python parallel_hycom_redownload_resilient.py merge --workers 4

Legacy CEC entrypoints remain supported. Pass --module-dir to select an
alternative hycom_download.py repair implementation.
"""

from __future__ import annotations

import argparse
import contextlib
import hashlib
import html
import os
import random
import shutil
import tempfile
import time
from pathlib import Path
from typing import Any, Callable, Iterator
from urllib.parse import parse_qsl, urlencode, urlsplit, urlunsplit

import numpy as np
import pandas as pd
import requests
import xarray as xr


# Manifest handling, worker recovery and CLI options share one implementation.
from .hycom_redownload import (
    ROW_ID_COLUMN,
    _build_parser,
    _load_repair_module,
    _run_worker,
    atomic_csv_write,
    ensure_beneath,
    exclusive_target_lock,
    merge,
    prepare,
    resolve_output_path,
    worker_manifest_path,
)


DEFAULT_HYCOM_DIR = Path("/Volumes/Compute/CEC/Data/hycom_currents_humboldt")
DEFAULT_DATA_ROOT = Path("/Volumes/Compute/CEC/Data")


@contextlib.contextmanager
def shared_download_slot(
    hycom_dir: Path,
    max_active: int,
    poll_seconds: float = 2.0,
    poll_jitter: float = 2.0,
) -> Iterator[int]:
    """Limit active HYCOM requests across all terminal processes.

    Each process tries to acquire one of ``max_active`` advisory lock files.
    The lock is held for an entire file replacement, including its retries and
    single-time fallback, so no more than ``max_active`` files are requesting
    HYCOM data simultaneously.
    """
    if max_active < 1:
        raise ValueError("--max-active-downloads must be at least 1")

    try:
        import fcntl
    except ImportError as exc:  # pragma: no cover - macOS/Linux only
        raise RuntimeError(
            "Shared download throttling requires the Unix fcntl module."
        ) from exc

    slot_dir = hycom_dir / ".hycom_download_slots"
    slot_dir.mkdir(parents=True, exist_ok=True)

    while True:
        indices = list(range(max_active))
        random.shuffle(indices)
        for index in indices:
            path = slot_dir / f"slot_{index}.lock"
            handle = open(path, "a+", encoding="utf-8")
            try:
                fcntl.flock(
                    handle.fileno(),
                    fcntl.LOCK_EX | fcntl.LOCK_NB,
                )
            except BlockingIOError:
                handle.close()
                continue

            try:
                handle.seek(0)
                handle.truncate()
                handle.write(
                    f"pid={os.getpid()} acquired={pd.Timestamp.now(tz='UTC').isoformat()}\n"
                )
                handle.flush()
                os.fsync(handle.fileno())
                yield index
            finally:
                with contextlib.suppress(OSError):
                    fcntl.flock(handle.fileno(), fcntl.LOCK_UN)
                handle.close()
            return

        delay = max(0.1, float(poll_seconds))
        delay += random.uniform(0.0, max(0.0, float(poll_jitter)))
        time.sleep(delay)


def import_repair_functions(module_dir: Path):
    """Return selected repair helpers, bypassing Python's module cache."""
    module = _load_repair_module(module_dir)
    return module.repair_missing_files, module._read_uv


def _manifest_target_rows(
    manifest_path: Path,
    target: Path,
    data_root: Path,
) -> tuple[pd.DataFrame, list[int]]:
    frame = pd.read_csv(manifest_path, dtype=str, keep_default_na=False)
    rows = [
        index
        for index, raw in enumerate(frame["output_path"])
        if resolve_output_path(raw, data_root) == target
    ]
    return frame, rows


def _find_recorded_ncss_url(row: pd.Series) -> str:
    """Return a complete bounded NCSS URL stored in a manifest row."""
    required = {"north", "south", "east", "west", "var"}
    for key in (
        "repair_url",
        "request_url",
        "download_url",
        "url",
        "source_url",
    ):
        raw = html.unescape(str(row.get(key, "")).strip())
        if not raw:
            continue
        parsed = urlsplit(raw)
        query = dict(parse_qsl(parsed.query, keep_blank_values=True))
        if (
            parsed.scheme in {"http", "https"}
            and "/ncss/" in parsed.path
            and required.issubset(query)
        ):
            return raw
    raise ValueError(
        "The manifest row does not contain a complete bounded NCSS request URL."
    )


def _single_time_url(base_url: str, timestamp: pd.Timestamp, accept: str) -> str:
    """Convert a time-range NCSS URL to one exact-time request."""
    parsed = urlsplit(base_url)
    remove = {
        "time",
        "time_start",
        "time_end",
        "time_duration",
        "timeStride",
        "temporal",
        "accept",
    }
    pairs = [
        (key, value)
        for key, value in parse_qsl(parsed.query, keep_blank_values=True)
        if key not in remove
    ]
    stamp = pd.Timestamp(timestamp)
    if stamp.tzinfo is not None:
        stamp = stamp.tz_convert("UTC").tz_localize(None)
    pairs.extend(
        [
            ("time", stamp.strftime("%Y-%m-%dT%H:%M:%SZ")),
            ("accept", accept),
        ]
    )
    return urlunsplit(
        (
            parsed.scheme,
            parsed.netloc,
            parsed.path,
            urlencode(pairs, doseq=True),
            "",
        )
    )


def _align_piece_to_template(
    frame: xr.Dataset,
    template: xr.Dataset,
    timestamp: pd.Timestamp,
) -> xr.Dataset:
    """Require the downloaded single-time grid to match the old file exactly."""
    actual_times = pd.DatetimeIndex(frame.time.values)
    expected = pd.Timestamp(timestamp)
    if expected.tzinfo is not None:
        expected = expected.tz_convert("UTC").tz_localize(None)

    if len(actual_times) == 0:
        raise ValueError("Single-time response contains no time coordinate.")
    deltas = np.abs(actual_times.values - expected.to_datetime64())
    index = int(np.argmin(deltas))
    if pd.Timedelta(deltas[index]) > pd.Timedelta(seconds=1):
        raise ValueError(
            f"Single-time response returned {actual_times[index]} instead of {expected}."
        )
    frame = frame.isel(time=slice(index, index + 1))

    for name in ("depth", "lat", "lon"):
        source = pd.Index(np.asarray(frame[name].values, dtype=float))
        wanted = np.asarray(template[name].values, dtype=float)
        indices = source.get_indexer(
            wanted,
            method="nearest",
            tolerance=1.0e-5,
        )
        if (indices < 0).any() or len(np.unique(indices)) != len(indices):
            raise ValueError(
                f"Single-time response {name} grid does not match the template."
            )
        frame = frame.isel({name: indices}).assign_coords({name: template[name].copy()})

    frame = frame.assign_coords(time=np.asarray([expected.to_datetime64()])).transpose(
        "time", "depth", "lat", "lon"
    )

    expected_shape = (
        1,
        template.sizes["depth"],
        template.sizes["lat"],
        template.sizes["lon"],
    )
    for variable in ("water_u", "water_v"):
        if tuple(frame[variable].shape) != expected_shape:
            raise ValueError(
                f"{variable} has shape {frame[variable].shape}; expected {expected_shape}."
            )
    if not (
        np.isfinite(frame.water_u.values) & np.isfinite(frame.water_v.values)
    ).any():
        raise ValueError("Single-time response has no finite paired U/V values.")
    return frame


def _download_single_time_piece(
    session: requests.Session,
    base_url: str,
    timestamp: pd.Timestamp,
    template: xr.Dataset,
    read_uv: Callable[[Path], xr.Dataset],
    target_dir: Path,
    *,
    retries: int,
    timeout: tuple[float, float],
    backoff_seconds: float,
) -> tuple[xr.Dataset, str]:
    """Download one native model time, preferring NetCDF-3 over NetCDF-4."""
    errors: list[str] = []
    for accept in ("netCDF", "netCDF4"):
        url = _single_time_url(base_url, timestamp, accept)
        for attempt in range(1, retries + 1):
            fd, temporary_name = tempfile.mkstemp(
                prefix="hycom_single_time.",
                suffix=".part.nc",
                dir=target_dir,
            )
            os.close(fd)
            temporary = Path(temporary_name)
            permanent_error = False
            retry_after = 0.0
            try:
                print(
                    "  single-time %s using %s (attempt %d/%d)"
                    % (
                        pd.Timestamp(timestamp).strftime("%Y-%m-%d %H:%MZ"),
                        accept,
                        attempt,
                        retries,
                    ),
                    flush=True,
                )
                with session.get(url, stream=True, timeout=timeout) as response:
                    permanent_error = (
                        400 <= response.status_code < 500
                        and response.status_code not in (408, 429)
                    )
                    try:
                        retry_after = max(
                            0.0,
                            float(response.headers.get("Retry-After", 0.0)),
                        )
                    except (TypeError, ValueError):
                        pass
                    response.raise_for_status()
                    with open(temporary, "wb") as output:
                        for chunk in response.iter_content(1024 * 1024):
                            if chunk:
                                output.write(chunk)
                        output.flush()
                        os.fsync(output.fileno())

                frame = read_uv(temporary)
                frame = _align_piece_to_template(
                    frame,
                    template,
                    timestamp,
                )
                return frame, accept
            except (
                requests.RequestException,
                OSError,
                ValueError,
                KeyError,
                IndexError,
            ) as exc:
                errors.append(f"{url}: {type(exc).__name__}: {exc}")
                if permanent_error:
                    break
                if attempt < retries:
                    delay = max(
                        retry_after,
                        float(backoff_seconds) * 2 ** (attempt - 1),
                    )
                    delay += random.uniform(0.0, max(1.0, delay * 0.25))
                    time.sleep(delay)
            finally:
                temporary.unlink(missing_ok=True)

    raise RuntimeError(
        "No valid single-time NCSS response could be downloaded for "
        f"{pd.Timestamp(timestamp)}.\n" + "\n".join(errors)
    )


def _commit_split_success(
    manifest_path: Path,
    target: Path,
    data_root: Path,
    *,
    source_url: str,
    template_path: Path,
    request_count: int,
    accepts: list[str],
    file_size: int,
    sha256: str,
) -> None:
    fresh, rows = _manifest_target_rows(manifest_path, target, data_root)
    if len(rows) != 1:
        raise ValueError(f"Expected one manifest row for {target}; found {len(rows)}.")
    index = rows[0]
    values: dict[str, Any] = {
        "status": "downloaded",
        "repair_status": "repaired",
        "repair_mode": "ncss_single_time",
        "repair_url": source_url,
        "repaired_at_utc": pd.Timestamp.now(tz="UTC").isoformat(),
        "file_size_bytes": str(file_size),
        "sha256": sha256,
        "repair_template": str(template_path),
        "repair_time_request_count": str(request_count),
        "repair_accept_formats": ",".join(sorted(set(accepts))),
        "repair_error": "",
    }
    for key, value in values.items():
        if key not in fresh.columns:
            fresh[key] = ""
        fresh.at[index, key] = value

    backup = manifest_path.with_suffix(manifest_path.suffix + ".bak")
    if not backup.exists():
        shutil.copy2(manifest_path, backup)
    atomic_csv_write(fresh, manifest_path)


def repair_by_single_times(
    target: Path,
    manifest_path: Path,
    data_root: Path,
    template_path: Path | None,
    read_uv: Callable[[Path], xr.Dataset],
    *,
    retries: int,
    timeout: tuple[float, float],
    backoff_seconds: float,
) -> None:
    """Replace one daily file with eight small NCSS time-slice requests."""
    if template_path is None or not template_path.is_file():
        raise ValueError(
            "Single-time fallback requires the original nonempty file as a template."
        )

    manifest, rows = _manifest_target_rows(manifest_path, target, data_root)
    if len(rows) != 1:
        raise ValueError(f"Expected one manifest row for {target}; found {len(rows)}.")
    row = manifest.iloc[rows[0]]
    base_url = _find_recorded_ncss_url(row)
    template = read_uv(template_path)
    expected_times = pd.DatetimeIndex(template.time.values)
    if expected_times.hasnans or expected_times.has_duplicates:
        raise ValueError("Template contains missing or duplicate timestamps.")
    if len(expected_times) == 0:
        raise ValueError("Template contains no timestamps.")

    pieces: list[xr.Dataset] = []
    accepts: list[str] = []
    with requests.Session() as session:
        session.headers.update({"User-Agent": "CEC-HYCOM-resilient-single-time/2.0"})
        for timestamp in expected_times:
            piece, accept = _download_single_time_piece(
                session,
                base_url,
                timestamp,
                template,
                read_uv,
                target.parent,
                retries=retries,
                timeout=timeout,
                backoff_seconds=backoff_seconds,
            )
            pieces.append(piece)
            accepts.append(accept)

    combined = xr.concat(
        pieces,
        dim="time",
        data_vars="all",
        coords="minimal",
        compat="override",
        join="exact",
    ).transpose("time", "depth", "lat", "lon")
    combined = combined.assign_coords(
        time=template.time.copy(),
        depth=template.depth.copy(),
        lat=template.lat.copy(),
        lon=template.lon.copy(),
    )
    combined.attrs.update(
        repair_source_url=base_url,
        repair_mode="ncss_single_time",
        repair_template=str(template_path),
    )
    for name in combined.variables:
        combined[name].encoding = {}

    fd, temporary_name = tempfile.mkstemp(
        prefix=target.stem + ".",
        suffix=".assembled.nc",
        dir=target.parent,
    )
    os.close(fd)
    temporary = Path(temporary_name)
    try:
        combined.to_netcdf(temporary)
        validated = read_uv(temporary)
        if not pd.DatetimeIndex(validated.time.values).equals(expected_times):
            raise ValueError("Assembled file timestamps differ from the template.")
        for name in ("depth", "lat", "lon"):
            if not np.allclose(
                validated[name].values,
                template[name].values,
                rtol=0.0,
                atol=1.0e-8,
                equal_nan=True,
            ):
                raise ValueError(
                    f"Assembled file {name} coordinate differs from the template."
                )
        if validated.water_u.shape != template.water_u.shape:
            raise ValueError("Assembled water_u shape differs from the template.")
        if validated.water_v.shape != template.water_v.shape:
            raise ValueError("Assembled water_v shape differs from the template.")
        if not (
            np.isfinite(validated.water_u.values)
            & np.isfinite(validated.water_v.values)
        ).any():
            raise ValueError("Assembled file has no finite paired U/V velocities.")

        digest = hashlib.sha256()
        with open(temporary, "rb") as handle:
            for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(chunk)
        file_size = temporary.stat().st_size
        os.replace(temporary, target)
    finally:
        temporary.unlink(missing_ok=True)

    _commit_split_success(
        manifest_path,
        target,
        data_root,
        source_url=base_url,
        template_path=template_path,
        request_count=len(pieces),
        accepts=accepts,
        file_size=file_size,
        sha256=digest.hexdigest(),
    )
    print(
        f"Repaired {target.name} with {len(pieces)} single-time requests; "
        "manifest updated.",
        flush=True,
    )


def repair_with_strategy(
    target: Path,
    manifest_path: Path,
    data_root: Path,
    template_path: Path | None,
    repair_missing_files: Callable[..., Any],
    read_uv: Callable[[Path], xr.Dataset],
    args: argparse.Namespace,
) -> str:
    """Use whole-day NCSS when suitable, with a single-time fallback."""
    requested_strategy = args.strategy
    strategy = requested_strategy
    if strategy == "auto":
        # HYCOM's old expt_53.X yearly aggregation is the source most prone to
        # long server-side requests. Start with native single times there.
        can_split = template_path is not None and template_path.is_file()
        strategy = (
            "single-time"
            if "gofs31r_expt53x" in target.name.lower() and can_split
            else "whole-day"
        )

    timeout = (args.connect_timeout, args.read_timeout)
    if strategy == "single-time":
        repair_by_single_times(
            target,
            manifest_path,
            data_root,
            template_path,
            read_uv,
            retries=args.split_retries,
            timeout=timeout,
            backoff_seconds=args.backoff,
        )
        return "single-time"

    try:
        repair_missing_files(
            [target],
            manifest_path=manifest_path,
            data_root=data_root,
            template_path=template_path,
            retries=args.retries,
            timeout=timeout,
            backoff_seconds=args.backoff,
        )
        return "whole-day"
    except Exception as primary_error:
        if requested_strategy == "whole-day" or not args.split_fallback:
            raise
        print(
            "Whole-day request failed; changing to native single-time requests.\n"
            f"Primary error: {primary_error}",
            flush=True,
        )
        try:
            repair_by_single_times(
                target,
                manifest_path,
                data_root,
                template_path,
                read_uv,
                retries=args.split_retries,
                timeout=timeout,
                backoff_seconds=args.backoff,
            )
            return "single-time-fallback"
        except Exception as fallback_error:
            raise RuntimeError(
                "Both the whole-day request and the single-time fallback failed.\n"
                f"Whole-day error: {primary_error}\n"
                f"Single-time error: {fallback_error}"
            ) from fallback_error


def worker(args: argparse.Namespace) -> int:
    hycom_dir = args.hycom_dir.expanduser().resolve()

    def repair_factory():
        repair_missing_files, read_uv = import_repair_functions(args.module_dir)

        def repair_target(target, manifest, data_root, template):
            if args.start_jitter > 0:
                time.sleep(random.uniform(0.0, args.start_jitter))
            with shared_download_slot(
                hycom_dir,
                args.max_active_downloads,
                poll_seconds=args.slot_poll,
                poll_jitter=args.slot_jitter,
            ) as slot:
                print(
                    f"Using shared HYCOM download slot {slot + 1}/"
                    f"{args.max_active_downloads}.",
                    flush=True,
                )
                return repair_with_strategy(
                    target,
                    manifest,
                    data_root,
                    template,
                    repair_missing_files,
                    read_uv,
                    args,
                )

        return repair_target

    def cooldown():
        if args.failure_cooldown > 0:
            delay = args.failure_cooldown + random.uniform(
                0.0, args.failure_cooldown * 0.5
            )
            print(f"Cooling down for {delay:.0f} seconds after failure.", flush=True)
            time.sleep(delay)

    return _run_worker(args, repair_factory, on_failure=cooldown)


def build_parser() -> argparse.ArgumentParser:
    parser, worker_parser = _build_parser(
        worker_func=worker,
        description=(
            "Split a HYCOM manifest among independent terminal workers, "
            "redownload each worker's files with shared throttling and single-time fallback, and merge safely."
        ),
        hycom_dir=DEFAULT_HYCOM_DIR,
        data_root=DEFAULT_DATA_ROOT,
        retries=1,
        connect_timeout=30,
        read_timeout=900,
        backoff=20,
        worker_help={
            "--retries": "Whole-day request attempts before fallback (default: 1).",
            "--read-timeout": "Socket read timeout in seconds (default: 900).",
            "--backoff": "Initial exponential retry backoff in seconds (default: 20).",
        },
    )
    worker_parser.add_argument(
        "--strategy",
        choices=("auto", "whole-day", "single-time"),
        default="auto",
        help=(
            "Download strategy. auto uses single-time requests first for "
            "GOFS expt_53.X and whole-day requests elsewhere (default: auto)."
        ),
    )
    worker_parser.add_argument(
        "--split-retries",
        type=int,
        default=4,
        help="Attempts for each native single-time request (default: 4).",
    )
    worker_parser.add_argument(
        "--max-active-downloads",
        type=int,
        default=3,
        help=(
            "Maximum simultaneous HYCOM file downloads across all worker "
            "terminals using this output directory (default: 3)."
        ),
    )
    worker_parser.add_argument("--slot-poll", type=float, default=2)
    worker_parser.add_argument("--slot-jitter", type=float, default=3)
    worker_parser.add_argument(
        "--start-jitter",
        type=float,
        default=15,
        help="Random delay before each file request, in seconds (default: 15).",
    )
    worker_parser.add_argument(
        "--failure-cooldown",
        type=float,
        default=60,
        help="Base cooldown after a failed file, in seconds (default: 60).",
    )
    worker_parser.add_argument(
        "--split-fallback",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Fall back to single-time requests after a whole-day failure.",
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
