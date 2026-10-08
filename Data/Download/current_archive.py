"""Download bounded current chunks concurrently into resumable CF NetCDF archives.

Dates are UTC and intervals are half-open: start <= time < end. Each archive
has a fixed request.json and a download_manifest.csv. Rerun the same command
to resume, or run --worker-index 0..N-1 in N terminals with --workers N.
All workers share file locks and atomically update the master manifest.

From the projects directory, for example:
    python -m HyperNav.Utilities.Data.current_archive copernicus \
        --start 2023-12-16 --end 2024-05-15 --workers 4 \
        --dataset-id cmems_mod_glo_phy_my_0.083deg_P1D-m --output-dir /data/copernicus
    python -m HyperNav.Utilities.Data.current_archive wcofs \
        --start YYYY-MM-DD --end YYYY-MM-DD --workers 4 --output-dir /data/wcofs

Authenticate Copernicus separately with `copernicusmarine login`, or its
documented environment variables. WCOFS defaults to the historical NCEI ROMS
archive with physical regridding. Use --source regulargrid for NOAA's short
rolling regular-grid archive.

Saved files contain U/V (and W only if actually supplied), in m/s, on
time/depth/lat/lon with physical depth in meters positive down. They can be
passed directly to CEC.current_load_year.ModelCurrent(files_by_year). Include
data bracketing the simulation dates, not just samples inside the season.
"""

from concurrent.futures import ThreadPoolExecutor, as_completed
from contextlib import contextmanager
import argparse
import hashlib
import json
from numbers import Integral
import os
from pathlib import Path
import re
import tempfile
import threading

import numpy as np
import pandas as pd
import xarray as xr


# HDF5/netCDF C libraries may be unsafe when called from simultaneous threads.
# Fetching and computing run concurrently; local NetCDF operations use this lock.
NETCDF_LOCK = threading.RLock()


@contextmanager
def _file_lock(path):
    import fcntl

    with path.open("a+") as handle:
        fcntl.flock(handle.fileno(), fcntl.LOCK_EX)
        try:
            yield
        finally:
            fcntl.flock(handle.fileno(), fcntl.LOCK_UN)


def _atomic_text(path, content):
    descriptor, temporary = tempfile.mkstemp(prefix=path.name + ".", dir=path.parent)
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8", newline="") as handle:
            handle.write(content)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
    finally:
        Path(temporary).unlink(missing_ok=True)


def _utc(value):
    stamp = pd.Timestamp(value)
    if pd.isna(stamp):
        raise ValueError("Use finite UTC dates or timestamps.")
    if stamp.tzinfo is not None:
        stamp = stamp.tz_convert("UTC").tz_localize(None)
    return stamp


def _validate_currents(dataset, start, end):
    if not isinstance(dataset, xr.Dataset):
        raise TypeError("Current fetcher must return an xarray.Dataset.")
    if not {"U", "V"}.issubset(dataset.data_vars):
        raise ValueError("Current chunks require U and V velocity variables.")
    dimensions = ("time", "depth", "lat", "lon")
    for name in dimensions:
        if (
            name not in dataset.coords
            or dataset[name].dims != (name,)
            or not dataset[name].size
        ):
            raise ValueError(f"{name} must be a nonempty one-dimensional coordinate.")
        values = dataset[name].values
        if name == "time":
            if not np.issubdtype(values.dtype, np.datetime64) or np.isnat(values).any():
                raise ValueError(
                    "Time must be decoded, nonmissing Gregorian datetime64."
                )
        elif (
            not np.issubdtype(values.dtype, np.number) or not np.isfinite(values).all()
        ):
            raise ValueError(f"{name} coordinates must be finite numbers.")
        if len(values) > 1 and not np.all(values[1:] > values[:-1]):
            raise ValueError(f"{name} coordinates must be unique and increasing.")
    times = dataset.time.values
    if times[0] < start.to_datetime64() or times[-1] >= end.to_datetime64():
        raise ValueError(
            "Current timestamps must lie within the requested half-open chunk."
        )
    if (
        dataset.depth.attrs.get("units") != "m"
        or dataset.depth.attrs.get("positive") != "down"
        or (dataset.depth < 0).any()
    ):
        raise ValueError(
            "Depth must be physical meters, nonnegative and positive down."
        )
    for name, units, limit in (
        ("lat", "degrees_north", 90),
        ("lon", "degrees_east", 180),
    ):
        if (
            dataset[name].attrs.get("units") != units
            or (np.abs(dataset[name]) > limit).any()
        ):
            raise ValueError(f"{name} must use geographic {units} coordinates.")
    if (np.diff(dataset.lon.values) > 180).any():
        raise ValueError("Regional date-line grids need regridding before archiving.")
    for name in ("U", "V", "W"):
        if name not in dataset:
            continue
        velocity = dataset[name]
        if velocity.dims != dimensions or velocity.attrs.get("units") not in (
            "m/s",
            "m s-1",
        ):
            raise ValueError(f"{name} must use m/s on the time/depth/lat/lon grid.")
        if (
            not np.issubdtype(velocity.dtype, np.number)
            or np.isinf(velocity.values).any()
        ):
            raise ValueError(
                f"{name} must contain finite physical velocities or NaN masks."
            )
        if name == "W" and (
            velocity.attrs.get("positive", "down") != "down"
            or velocity.attrs.get("standard_name") == "upward_sea_water_velocity"
        ):
            raise ValueError("W must be positive downward.")
    if not (np.isfinite(dataset.U.values) & np.isfinite(dataset.V.values)).any():
        raise ValueError("Current chunk contains no jointly valid U/V samples.")
    return {
        "records": len(times),
        "first_time": pd.Timestamp(times[0]).isoformat(),
        "last_time": pd.Timestamp(times[-1]).isoformat(),
    }


def download_record(
    fetch_chunk,
    *,
    start,
    end,
    output_dir,
    model,
    workers=4,
    chunk_days=1,
    retries=3,
    worker_index=None,
    settings=None,
    engine=None,
):
    """Fetch canonical current chunks, preserving successes and recording failures.

    fetch_chunk(start, end) receives timezone-naive UTC pandas Timestamps.
    retries is the maximum number of attempts per chunk. worker_index selects
    every Nth chunk for one terminal/process; otherwise workers is the thread
    count. Return the complete master manifest, including other workers' rows.
    Local file locks make duplicated/restarted workers safe on macOS/Linux.
    """
    start, end = _utc(start), _utc(end)
    if end <= start:
        raise ValueError("end must be later than start.")
    for name, value in (
        ("workers", workers),
        ("chunk_days", chunk_days),
        ("retries", retries),
    ):
        if isinstance(value, bool) or not isinstance(value, Integral) or value < 1:
            raise ValueError(f"{name} must be a positive integer.")
    if worker_index is not None and (
        isinstance(worker_index, bool)
        or not isinstance(worker_index, Integral)
        or not 0 <= worker_index < workers
    ):
        raise ValueError("worker_index must lie between zero and workers - 1.")
    if not callable(fetch_chunk):
        raise TypeError("fetch_chunk must be callable.")
    if not isinstance(model, str) or not re.fullmatch(
        r"[A-Za-z0-9][A-Za-z0-9_.-]{0,79}", model
    ):
        raise ValueError("model must be a simple filename-safe label.")
    configuration = {
        "model": model,
        "start": start.isoformat(),
        "end": end.isoformat(),
        "chunk_days": int(chunk_days),
        "settings": dict(settings or {}),
        "engine": engine,
    }
    encoded = json.dumps(configuration, sort_keys=True, allow_nan=False)
    configuration = json.loads(encoded)
    fingerprint = hashlib.sha256(encoded.encode()).hexdigest()
    directory = Path(output_dir).expanduser().resolve()
    directory.mkdir(parents=True, exist_ok=True)
    manifest_path = directory / "download_manifest.csv"
    lock_path = directory / "download_manifest.lock"
    request_path = directory / "request.json"
    rows = []
    chunk_start = start
    while chunk_start < end:
        chunk_end = min(chunk_start + pd.Timedelta(days=int(chunk_days)), end)
        index = len(rows)
        filename = f"{model}_{index:05d}_{chunk_start:%Y%m%dT%H%M%S}_{chunk_end:%Y%m%dT%H%M%S}.nc"
        rows.append(
            {
                "chunk_id": index,
                "start": chunk_start.isoformat(),
                "end": chunk_end.isoformat(),
                "output_path": str(directory / filename),
                "status": "pending",
                "error": "",
                "attempts": 0,
                "records": 0,
                "first_time": "",
                "last_time": "",
            }
        )
        chunk_start = chunk_end
    expected = pd.DataFrame(rows)
    with _file_lock(lock_path):
        if request_path.exists():
            if json.loads(request_path.read_text()) != configuration:
                raise ValueError(
                    "Archive request differs from this directory's saved settings; use a new output directory."
                )
        elif manifest_path.exists() or any(directory.glob("*.nc")):
            raise ValueError(
                "Existing archive lacks request.json; use a new output directory."
            )
        else:
            _atomic_text(request_path, json.dumps(configuration, indent=2) + "\n")
        if manifest_path.exists():
            saved = pd.read_csv(manifest_path, keep_default_na=False)
            keys = ["chunk_id", "start", "end", "output_path"]
            if not set(expected.columns).issubset(saved.columns) or not saved[
                keys
            ].equals(expected[keys]):
                raise ValueError(
                    "Archive manifest does not match the saved chunk schedule."
                )
        else:
            _atomic_text(manifest_path, expected.to_csv(index=False))

    def update(index, values):
        with _file_lock(lock_path):
            saved = pd.read_csv(manifest_path, keep_default_na=False)
            for key, value in values.items():
                saved.loc[saved.chunk_id == index, key] = value
            _atomic_text(manifest_path, saved.to_csv(index=False))

    def run_chunk(row):
        index = row["chunk_id"]
        left, right = _utc(row["start"]), _utc(row["end"])
        target = Path(row["output_path"])
        with _file_lock(target.with_suffix(".nc.lock")):
            if target.is_file():
                try:
                    with NETCDF_LOCK, xr.open_dataset(target) as cached:
                        if cached.attrs.get("archive_settings_sha256") != fingerprint:
                            raise ValueError(
                                "Cached chunk has different archive settings."
                            )
                        if (
                            cached.attrs.get("requested_start") != left.isoformat()
                            or cached.attrs.get("requested_end") != right.isoformat()
                        ):
                            raise ValueError(
                                "Cached chunk has different requested dates."
                            )
                        metrics = _validate_currents(cached, left, right)
                    update(index, {"status": "downloaded", "error": "", **metrics})
                    return
                except Exception:
                    # Refetch through a separate temporary file; preserve the old
                    # bytes if the request or replacement validation fails.
                    pass
            error = ""
            for attempt in range(1, int(retries) + 1):
                dataset = None
                source = None
                temporary = None
                try:
                    source = dataset = fetch_chunk(left, right)
                    if not isinstance(dataset, xr.Dataset):
                        raise TypeError(
                            "Current fetcher must return an xarray.Dataset."
                        )
                    dataset.load()
                    metrics = _validate_currents(dataset, left, right)
                    dataset = dataset[
                        [name for name in ("U", "V", "W") if name in dataset]
                    ].copy()
                    dataset.attrs.update(
                        archive_settings_sha256=fingerprint,
                        requested_start=left.isoformat(),
                        requested_end=right.isoformat(),
                        model=model,
                        Conventions="CF-1.8",
                    )
                    for variable in dataset.variables.values():
                        variable.encoding.clear()
                    descriptor, temporary_name = tempfile.mkstemp(
                        prefix=target.name + ".", suffix=".nc", dir=directory
                    )
                    os.close(descriptor)
                    temporary = Path(temporary_name)
                    with NETCDF_LOCK:
                        dataset.to_netcdf(temporary, engine=engine)
                        with xr.open_dataset(temporary, engine=engine) as written:
                            _validate_currents(written, left, right)
                    os.replace(temporary, target)
                    update(
                        index,
                        {
                            "status": "downloaded",
                            "error": "",
                            "attempts": attempt,
                            **metrics,
                        },
                    )
                    print(f"{target.name}: downloaded", flush=True)
                    return
                except Exception as exc:
                    error = f"{type(exc).__name__}: {exc}"
                finally:
                    if dataset is not None and hasattr(dataset, "close"):
                        dataset.close()
                    if (
                        source is not None
                        and source is not dataset
                        and hasattr(source, "close")
                    ):
                        source.close()
                    if temporary is not None:
                        temporary.unlink(missing_ok=True)
            update(
                index,
                {
                    "status": "failed",
                    "error": error,
                    "attempts": int(retries),
                    "records": 0,
                    "first_time": "",
                    "last_time": "",
                },
            )
            print(f"{target.name}: failed ({error})", flush=True)

    assigned = rows if worker_index is None else rows[int(worker_index) :: int(workers)]
    with ThreadPoolExecutor(
        max_workers=int(workers) if worker_index is None else 1
    ) as pool:
        futures = [pool.submit(run_chunk, row) for row in assigned]
        for future in as_completed(futures):
            future.result()
    with _file_lock(lock_path):
        return pd.read_csv(manifest_path, keep_default_na=False)
