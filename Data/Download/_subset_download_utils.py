"""Small I/O helpers shared by regional model subset downloaders."""

from pathlib import Path
import os
import time
import xml.etree.ElementTree as ET
import requests
from tempfile import TemporaryDirectory, mkstemp

import xarray as xr

NETCDF_SIGNATURES = (b"\x89HDF\r\n\x1a\n", b"CDF\x01", b"CDF\x02", b"CDF\x05")


def read_netcdf_bytes(content, *, prefix="model_subset_"):
    """Decode a complete NetCDF3/4 response and detach it from temporary files."""
    if not content.startswith(NETCDF_SIGNATURES):
        raise OSError(f"Not NetCDF: {content[:600]!r}")
    with TemporaryDirectory(prefix=prefix) as folder:
        path = Path(folder) / "subset.nc"
        path.write_bytes(content)
        with xr.open_dataset(path, engine="netcdf4") as source:
            return source.load()


def save_netcdf_atomic(dataset, path, *, encode_time=False):
    """Publish a NetCDF only after the entire staged write succeeds."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary = mkstemp(prefix=f".{path.name}.", suffix=".nc", dir=path.parent)
    os.close(descriptor)
    try:
        saved = dataset.copy(deep=False)
        if encode_time:
            # Isolate encoding changes from the caller's coordinate variable.
            saved["time"] = dataset.time.copy(deep=False)
            saved.time.encoding = {
                "dtype": "int64",
                "units": "hours since 2000-01-01 00:00:00",
                "calendar": dataset.time.encoding.get("calendar", "proleptic_gregorian"),
            }
        saved.to_netcdf(temporary, engine="netcdf4")
        os.replace(temporary, path)
    finally:
        Path(temporary).unlink(missing_ok=True)


def request_bytes(
    session, url, *, params=None, timeout=120, max_attempts=3, retry_wait=2,
    reader=None, max_bytes=None, connect_timeout=20, retry_cap=15,
    failure_message=None, retry_reporter=None,
):
    """Retry transport/decoding failures without retrying permanent HTTP 4xx.

    A bounded streamed response is used when ``max_bytes`` is provided; model
    engines retain their own connection timeout, retry cap, and decoder checks.
    """
    for attempt in range(max_attempts):
        try:
            options = {"params": params, "timeout": (connect_timeout, timeout)}
            if max_bytes is not None:
                options["stream"] = True
            with session.get(url, **options) as response:
                response.raise_for_status()
                length = response.headers.get("Content-Length")
                uncompressed_length = length and not response.headers.get("Content-Encoding")
                if max_bytes is None:
                    content = response.content
                else:
                    if uncompressed_length and int(length) > max_bytes:
                        raise ValueError("NCSS response exceeds the expected regional subset size.")
                    chunks, size = [], 0
                    for chunk in response.iter_content(chunk_size=65536):
                        size += len(chunk)
                        if size > max_bytes:
                            raise ValueError("NCSS response exceeds the expected regional subset size.")
                        chunks.append(chunk)
                    content = b"".join(chunks)
                if uncompressed_length and len(content) != int(length):
                    raise OSError("Incomplete NCSS HTTP response.")
            return content if reader is None else reader(content)
        except (requests.RequestException, OSError, ValueError, RuntimeError, ET.ParseError) as error:
            if isinstance(error, requests.HTTPError) and error.response is not None:
                status = error.response.status_code
                if 400 <= status < 500 and status not in {408, 425, 429}:
                    raise
            if attempt + 1 == max_attempts:
                message = failure_message or f"HYCOM request failed after {max_attempts} attempts: {url}"
                raise RuntimeError(message) from error
            delay = min(float(retry_wait) * 2**attempt, retry_cap)
            if retry_reporter is None:
                print(
                    f"  HYCOM attempt {attempt + 1}/{max_attempts} failed: "
                    f"{type(error).__name__}: {str(error)[:240]}. "
                    f"Retrying in {delay:g}s (read timeout {timeout:g}s).", flush=True,
                )
            else:
                retry_reporter(error, attempt, delay)
            time.sleep(delay)
