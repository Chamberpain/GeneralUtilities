"""Daily, resumable HYCOM October current subsets from verified archives.

1994–2015 use GOFS 3.1 reanalysis; October 2016–2023 use the documented
GOFS 3.1 analysis experiments; 2024 onward use ESPC-D-V02 annual analysis.
The public NCSS service limits requests to one day. Older archives combine
U/V; ESPC archives them separately. Every source timestamp and grid coordinate
is checked before publication. Missing ocean values remain NaN.
"""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import tempfile
import time
import xml.etree.ElementTree as ET

import numpy as np
import pandas as pd
import requests
import xarray as xr

from ._subset_download_utils import read_netcdf_bytes, save_netcdf_atomic, request_bytes
from .hycom_download import opendap_subset as _opendap_subset


BASE_URL = "https://ncss.hycom.org/thredds/ncss/grid/ESPC-D-V02"
PROCESSOR = "hycom-october-archive-v2"
COMPONENTS = {"u3z": "water_u", "v3z": "water_v"}


def archive_sources(year):
    """Describe catalogue-verified sources for October of a supported year.

    The selector describes experiment routing, not availability of every date.
    `_archive_axes` checks the actual timestamps before a current request.
    2016–2018 experiments cover only parts of their respective calendar years.
    GOFS's 41-layer model output is interpolated to 40 standard physical depths;
    the exact exported depth coordinates are always read from source metadata.
    """
    if isinstance(year, bool) or not isinstance(year, (int, np.integer)) or year < 1994:
        raise ValueError("HYCOM October year must be an integer from 1994 onward.")
    year = int(year)
    ncss = "https://ncss.hycom.org/thredds/ncss/grid"
    catalogue = "https://tds.hycom.org/thredds/catalogs"
    if year >= 2024:
        family = "ESPC-D-V02"
        model_name = "HYCOM ESPC-D-V02 annual analysis"
        paths = [f"ESPC-D-V02/{component}/{year}" for component in COMPONENTS]
        groups = [[variable] for variable in COMPONENTS.values()]
        catalogue_url = f"{catalogue}/ESPC-D-V02.html"
    else:
        if year <= 2015:
            family = "GOFS-3.1-reanalysis"
            model_name = "HYCOM GOFS 3.1 reanalysis"
            path = f"GLBv0.08/expt_53.X/data/{year}"
            catalogue_url = f"{catalogue}/GLBv0.08/expt_53.X.html"
        else:
            model_name = "HYCOM GOFS 3.1 analysis"
            if year == 2016:
                experiment = "GLBv0.08/expt_57.2"
                path = experiment
            elif year == 2017:
                experiment = "GLBv0.08/expt_92.9"
                path = f"{experiment}/uv3z"
            elif year == 2018:
                experiment = "GLBu0.08/expt_93.0"
                path = f"{experiment}/uv3z"
            else:
                experiment = "GLBy0.08/expt_93.0"
                path = f"{experiment}/uv3z/{year}"
            family = f"GOFS-3.1-{experiment.replace('/', '-')}"
            suffix = "_uv3z" if year >= 2019 else ""
            catalogue_url = f"{catalogue}/{experiment}{suffix}.html"
        paths = [path]
        groups = [list(COMPONENTS.values())]
        if year == 2017:
            # expt_57.7 ends on October 1 at 09Z; expt_92.9 starts at 12Z.
            # Both archives are needed for the first day's complete native grid.
            paths.insert(0, "GLBv0.08/expt_57.7")
            groups.insert(0, list(COMPONENTS.values()))
            family = "GOFS-3.1-GLBv0.08-expt_57.7+92.9"
            catalogue_url = f"{catalogue}/GLBv0.08/expt_57.7.html\n{catalogue_url}"
    return {
        "archive_family": family,
        "archive_year": year,
        "model_name": model_name,
        "source_urls": [f"{ncss}/{path}" for path in paths],
        "variable_groups": groups,
        "native_time_step_hours": 3,
        "catalogue_url": catalogue_url,
    }


# The driver uses the descriptive name for provenance; both names share one API.
describe_hycom_source = archive_sources


def _request(
    session,
    url,
    *,
    params=None,
    timeout=120,
    max_attempts=3,
    retry_wait=2,
    reader=None,
    max_bytes=2_000_000,
):
    """Retry transport, truncated responses, and invalid decoded responses."""
    return request_bytes(
        session, url, params=params, timeout=timeout, max_attempts=max_attempts,
        retry_wait=retry_wait, reader=reader, max_bytes=max_bytes,
    )


def _archive_axes(session, year, **request_options):
    """Read exact source axes; do not infer availability from a date label."""
    archive = archive_sources(year)
    reference = None
    source_times = {}
    source_axes = {}
    variable_times = {
        name: pd.DatetimeIndex([], tz="UTC") for name in COMPONENTS.values()
    }
    for url, variables in zip(
        archive["source_urls"], archive["variable_groups"], strict=True
    ):
        root = _request(
            session,
            f"{url}/dataset.xml",
            reader=ET.fromstring,
            **request_options,
        )
        axes = {}
        native_longitudes = None
        for axis in root.iter():
            if axis.tag.rsplit("}", 1)[-1] != "axis":
                continue
            name = axis.get("name")
            if name not in {"time", "depth", "lat", "lon"}:
                continue
            attributes = {
                node.get("name"): node.get("value")
                for node in axis
                if node.tag.rsplit("}", 1)[-1] == "attribute"
            }
            values = next(
                (node for node in axis if node.tag.rsplit("}", 1)[-1] == "values"), None
            )
            if values is None:
                raise ValueError(f"HYCOM metadata lacks {name} axis values.")
            if (values.text or "").strip():
                numbers = np.array((values.text or "").split(), dtype=float)
            else:
                count = int(values.get("npts", "0"))
                if not 0 < count <= 1_000_000:
                    raise ValueError(f"Invalid {name} metadata axis length.")
                numbers = float(values.get("start")) + float(
                    values.get("increment")
                ) * np.arange(count)
            if (
                numbers.size != int(axis.get("shape", "0"))
                or not np.isfinite(numbers).all()
            ):
                raise ValueError(f"Invalid {name} metadata axis values.")
            if name == "time":
                calendar = attributes.get("calendar", "standard")
                if calendar not in {"standard", "gregorian", "proleptic_gregorian"}:
                    raise ValueError(f"Unsupported HYCOM calendar {calendar!r}.")
                decoded = xr.coding.times.decode_cf_datetime(
                    numbers, attributes["units"], calendar, time_unit="ns"
                )
                axes[name] = pd.DatetimeIndex(decoded).tz_localize("UTC")
                if (
                    axes[name].hasnans
                    or not axes[name].is_unique
                    or not axes[name].is_monotonic_increasing
                ):
                    raise ValueError("HYCOM metadata times must be unique and sorted.")
            else:
                if not numbers.size or np.any(np.diff(numbers) <= 0):
                    raise ValueError(
                        f"HYCOM metadata {name} must be unique and increasing."
                    )
                if name == "depth" and (
                    attributes.get("units", "m") not in {"m", "meters", "metres"}
                    or attributes.get("positive", "down") != "down"
                ):
                    raise ValueError(
                        "HYCOM metadata depths must be positive-down metres."
                    )
                if name == "lon":
                    native_longitudes = numbers.copy()
                    # Consecutive GOFS experiments use -180..180 and 0..360
                    # conventions for the same physical nodes.
                    numbers = np.sort((numbers + 180) % 360 - 180)
                    if np.any(np.diff(numbers) <= 0):
                        raise ValueError(
                            "HYCOM geographic longitude nodes must be unique."
                        )
                axes[name] = numbers
        if set(axes) != {"time", "depth", "lat", "lon"}:
            raise ValueError("HYCOM metadata must provide time/depth/lat/lon axes.")
        if reference is not None:
            for name in ("depth", "lat", "lon"):
                if not np.array_equal(axes[name], reference[name]):
                    raise ValueError(f"HYCOM U/V {name} metadata differ.")
        for variable in variables:
            variable_times[variable] = variable_times[variable].union(axes["time"])
        source_times[url] = axes["time"]
        source_axes[url] = {**axes, "lon": native_longitudes}
        reference = axes
    if not variable_times["water_u"].equals(variable_times["water_v"]):
        raise ValueError("HYCOM U/V time metadata differ.")
    reference["time"] = variable_times["water_u"]
    reference["source_times"] = source_times
    reference["source_axes"] = source_axes
    return reference



def hycom_archive_times(year, timeout=120):
    """Return exact UTC archive timestamps for the year's selected experiment.

    Raise if separate U/V archives advertise different timestamps or grid axes.
    The selected experiment can cover only part of the year or extend past it.
    """
    if not np.isfinite(timeout) or timeout <= 0:
        raise ValueError("timeout must be finite and positive.")
    with requests.Session() as session:
        return _archive_axes(session, year, timeout=timeout)["time"]


def _validate(dataset, *, expected_times, expected_depths, expected_grid, variables):
    """Validate source/cached arrays and normalize geographic longitude."""
    required = {"time", "depth", "lat", "lon"}
    if not set(variables).issubset(dataset.data_vars) or any(
        set(dataset[name].dims) != required for name in variables
    ):
        raise ValueError("HYCOM fields must have time/depth/lat/lon dimensions.")
    if any(dataset[name].dims != (name,) for name in required):
        raise ValueError("HYCOM coordinates must be one-dimensional axes.")
    if not np.issubdtype(dataset.time.dtype, np.datetime64):
        raise ValueError("HYCOM times must decode to datetime64.")
    result = dataset[list(variables)].reset_coords(drop=True)
    result = result.assign_coords(
        lon=(result.lon.astype(float) + 180) % 360 - 180
    ).sortby("lon")
    actual = pd.DatetimeIndex(result.time.values).tz_localize("UTC")
    if not actual.equals(expected_times):
        raise ValueError(
            "HYCOM returned incomplete, duplicated, or unexpected timestamps."
        )
    for name, expected in {"depth": expected_depths, **expected_grid}.items():
        values = np.asarray(result[name], dtype=float)
        if values.shape != np.asarray(expected).shape or not np.allclose(
            values, expected, rtol=0, atol=1e-6
        ):
            raise ValueError(
                f"HYCOM {name} grid differs from the requested source grid."
            )
    if result.depth.attrs.get(
        "positive", "down"
    ).lower() != "down" or result.depth.attrs.get("units", "m").lower() not in {
        "m",
        "meter",
        "meters",
        "metre",
        "metres",
    }:
        raise ValueError("HYCOM depths must be positive-down metres.")
    for name in variables:
        units = result[name].attrs.get("units", "").lower().replace(" ", "")
        if units not in {"m/s", "ms-1", "ms^-1", "m.s-1"}:
            raise ValueError(f"HYCOM {name} must use metres per second.")
        if np.isinf(result[name].values).any():
            raise ValueError("HYCOM velocities contain infinity.")
    return result.transpose("time", "depth", "lat", "lon").load()


def _save(dataset, path):
    save_netcdf_atomic(dataset, path, encode_time=True)




def _retryable_transport_failure(error):
    """Distinguish server/transport failures from invalid scientific data."""
    seen = set()
    while error is not None and id(error) not in seen:
        seen.add(id(error))
        if isinstance(error, requests.HTTPError):
            return error.response is not None and (
                error.response.status_code in {408, 413, 425, 429}
                or error.response.status_code >= 500
            )
        if isinstance(error, requests.RequestException):
            return True
        if isinstance(error, OSError) and "Incomplete NCSS HTTP response" in str(error):
            return True
        error = error.__cause__ or error.__context__
    return False


def _fetch_fragment(
    session,
    url,
    *,
    variables,
    times,
    depths,
    retained,
    grid,
    cache,
    signature_data,
    signature,
    force,
    method,
    request_options,
    source_axes=None,
):
    """Cache a component packet immediately; shrink transport failures serially."""
    spec = {
        "url": url,
        "variables": variables,
        "times": [time.isoformat() for time in times],
    }
    fingerprint = hashlib.sha256(json.dumps(spec, sort_keys=True).encode()).hexdigest()
    path = cache / "packets" / f"{fingerprint[:24]}.nc"
    validation = {
        "expected_times": times,
        "expected_depths": retained,
        "expected_grid": grid,
        "variables": variables,
    }
    if path.is_file() and not force:
        try:
            with xr.open_dataset(path, engine="netcdf4") as source:
                if (
                    source.attrs.get("request_signature") != signature
                    or source.attrs.get("packet_fingerprint") != fingerprint
                    or any(
                        source.attrs.get(key) != value
                        for key, value in signature_data.items()
                    )
                ):
                    raise ValueError("HYCOM packet cache belongs to another request.")
                return _validate(source, **validation)
        except (OSError, ValueError, RuntimeError):
            pass
    label = f"{', '.join(variables)} {times[0]} → {times[-1]} ({len(times)} times)"
    split_path = path.with_suffix(".split.json")
    split_record = {
        "request_signature": signature,
        "packet_fingerprint": fingerprint,
        "time_count": len(times),
    }
    resume_split = False
    if len(times) > 1 and split_path.is_file() and not force:
        try:
            resume_split = json.loads(split_path.read_text()) == split_record
        except (OSError, ValueError):
            pass
    result, failure = None, None
    used_method = method
    if not resume_split and method in {"auto", "opendap"}:
        try:
            print(f"HYCOM OPeNDAP {label}", flush=True)
            opendap_options = (
                {} if source_axes is None else {"source_axes": source_axes}
            )
            result = _opendap_subset(
                session,
                url,
                variables=variables,
                times=times,
                depths=retained,
                grid=grid,
                timeout=request_options["timeout"],
                **opendap_options,
            )
            result = _validate(result, **validation)
            used_method = "opendap"
        except (
            ImportError,
            requests.RequestException,
            OSError,
            ValueError,
            RuntimeError,
        ) as error:
            failure = error
            if method == "auto":
                print(
                    f"  OPeNDAP failed ({type(error).__name__}: {str(error)[:200]}); trying NCSS.",
                    flush=True,
                )
    if result is None and not resume_split and method in {"auto", "ncss"}:

        def decode(content):
            decoded = _validate(
                read_netcdf_bytes(content, prefix="hycom_archive_"),
                **{**validation, "expected_depths": depths},
            )
            return decoded.sel(depth=retained)

        params = {
            "west": float(grid["lon"][0]),
            "east": float(grid["lon"][-1]),
            "south": float(grid["lat"][0]),
            "north": float(grid["lat"][-1]),
            "var": variables[0] if len(variables) == 1 else variables,
            "time_start": times[0].isoformat(),
            "time_end": times[-1].isoformat(),
            "timeStride": 1,
            "horizStride": 1,
            "vertStride": 1,
            "accept": "netCDF4",
            "addLatLon": "true",
        }
        maximum = int(
            len(times)
            * len(depths)
            * np.prod([len(v) for v in grid.values()])
            * len(variables)
            * 8
            + 2_000_000
        )
        try:
            print(
                f"HYCOM NCSS {label} (read timeout {request_options['timeout']:g}s)",
                flush=True,
            )
            result = _request(
                session,
                url,
                params=params,
                reader=decode,
                max_bytes=maximum,
                **request_options,
            )
            used_method = "ncss"
        except (requests.RequestException, OSError, ValueError, RuntimeError) as error:
            failure = error
    if result is None:
        if (
            resume_split
            or (failure is not None and _retryable_transport_failure(failure))
        ) and len(times) > 1:
            midpoint = len(times) // 2
            if not resume_split:
                split_path.parent.mkdir(parents=True, exist_ok=True)
                descriptor, temporary = tempfile.mkstemp(
                    prefix=".split_", suffix=".json", dir=split_path.parent
                )
                try:
                    with os.fdopen(descriptor, "w") as handle:
                        json.dump(split_record, handle)
                    os.replace(temporary, split_path)
                finally:
                    Path(temporary).unlink(missing_ok=True)
            print(f"  Retrying {label} as two smaller serial packets.", flush=True)
            parts = [
                _fetch_fragment(
                    session,
                    url,
                    variables=variables,
                    times=part,
                    depths=depths,
                    retained=retained,
                    grid=grid,
                    cache=cache,
                    signature_data=signature_data,
                    signature=signature,
                    force=force,
                    method=method,
                    request_options=request_options,
                    source_axes=source_axes,
                )
                for part in (times[:midpoint], times[midpoint:])
            ]
            result = xr.concat(
                parts,
                dim="time",
                data_vars="minimal",
                coords="minimal",
                compat="equals",
                join="exact",
                combine_attrs="drop_conflicts",
            )
            result = _validate(result, **validation)
            used_method = "recovered smaller packets"
        else:
            if failure is not None:
                failure.add_note(
                    f"HYCOM packet failed: {label}. Successful packet caches are retained in {cache}."
                )
                raise failure
            raise RuntimeError(f"No HYCOM data returned for {label}.")
    result.attrs.update(
        signature_data,
        request_signature=signature,
        packet_fingerprint=fingerprint,
        download_method=used_method,
    )
    _save(result, path)
    return result



def download_hycom_october(
    start,
    end,
    *,
    bounds,
    max_depth=500,
    output_path,
    cache_dir,
    force=False,
    timeout=120,
    max_attempts=3,
    retry_wait=2,
    chunk_hours=24,
    method="ncss",
):
    """Download a complete selected-archive interval, inclusive of both ends.

    Dates must fall in one calendar year and on the native three-hour grid.
    Requests are serial and cover at most chunk_hours (3..24) hours. NCSS
    fetches all layers; OPeNDAP fetches just the required native levels.
    Transport failures shrink to individual timestamps; invalid data raises.
    Successful component packets and daily caches are reusable after a failure.
    Final output replaces an existing file only
    after every expected timestamp, component, and grid has been verified.
    """
    start, end = pd.to_datetime(start, utc=True), pd.to_datetime(end, utc=True)
    if method not in {"auto", "ncss", "opendap"}:
        raise ValueError("HYCOM method must be auto, ncss, or opendap.")
    if (
        isinstance(chunk_hours, (bool, np.bool_))
        or not np.isscalar(chunk_hours)
        or not np.isfinite(chunk_hours)
        or not 3 <= chunk_hours <= 24
        or chunk_hours % 3
    ):
        raise ValueError("HYCOM chunk_hours must be a multiple of 3 between 3 and 24.")
    if pd.isna(start) or pd.isna(end) or end <= start or start.year != end.year:
        raise ValueError(
            "Dates must define an increasing interval within one calendar year."
        )
    if start != start.floor("3h") or end != end.floor("3h"):
        raise ValueError("HYCOM dates must align with the three-hour source grid.")
    if set(bounds) != {"west", "east", "south", "north"}:
        raise ValueError("bounds must contain west/east/south/north.")
    settings = {
        name: float(bounds[name]) for name in ("west", "east", "south", "north")
    }
    if not np.isfinite(
        [*settings.values(), max_depth, timeout, retry_wait]
    ).all() or not (
        -180 <= settings["west"] < settings["east"] <= 180
        and -90 <= settings["south"] < settings["north"] <= 90
    ):
        raise ValueError(
            "Provide finite geographic bounds without a dateline crossing."
        )
    if (
        max_depth < 0
        or timeout <= 0
        or retry_wait < 0
        or isinstance(max_attempts, bool)
        or not isinstance(max_attempts, int)
        or max_attempts < 1
    ):
        raise ValueError(
            "Require nonnegative depth/retry wait and positive timeout/attempt count."
        )
    expected_times = pd.date_range(start, end, freq="3h")
    archive = archive_sources(start.year)
    signature_data = {
        "processor": PROCESSOR,
        "archive_family": archive["archive_family"],
        "archive_year": archive["archive_year"],
        "model": archive["model_name"],
        "source_url": "\n".join(archive["source_urls"]),
        "source_catalogue": archive["catalogue_url"],
        "native_time_step_hours": archive["native_time_step_hours"],
        **settings,
        "max_depth": float(max_depth),
    }
    signature = hashlib.sha256(
        json.dumps(signature_data, sort_keys=True).encode()
    ).hexdigest()
    output_path, cache_dir = (
        Path(output_path).expanduser(),
        Path(cache_dir).expanduser(),
    )
    if output_path.is_dir() or cache_dir.is_file():
        raise ValueError(
            "output_path must name a file and cache_dir must name a directory."
        )
    cache = cache_dir / signature[:20]
    request_options = {
        "timeout": timeout,
        "max_attempts": max_attempts,
        "retry_wait": retry_wait,
    }
    with requests.Session() as session:
        session.headers.update({"User-Agent": "SeatrecPlots-HYCOM-Archive/1.0"})
        axes = _archive_axes(session, start.year, **request_options)
        missing_times = expected_times.difference(axes["time"])
        if not missing_times.empty:
            sample = ", ".join(timestamp.isoformat() for timestamp in missing_times[:8])
            extra = (
                f" (plus {len(missing_times) - 8} more)"
                if len(missing_times) > 8
                else ""
            )
            message = (
                "Requested HYCOM interval is not fully available in the analysis archive. "
                f"Missing {len(missing_times)} native 3-hour timestamp(s), first "
                f"{missing_times[0].isoformat()}, last {missing_times[-1].isoformat()}: {sample}{extra}. "
                f"Archive: {archive['model_name']}; available UTC range "
                f"{axes['time'][0].isoformat()} to {axes['time'][-1].isoformat()}. "
                "Select a gap-free interval; changing NCSS/OPeNDAP transport cannot fill missing source times."
            )
            if archive["model_name"] == "HYCOM GOFS 3.1 analysis":
                message += (
                    " Known GOFS analysis gaps: "
                    "https://www.hycom.org/faqs/477-did-you-know-there-are-missing-days-in-the-gofs-31-global-analysis"
                )
            raise ValueError(message)
        depths = axes["depth"]
        if depths[0] != 0 or max_depth > depths[-1]:
            raise ValueError("Requested HYCOM depth range is not available.")
        retained = depths[: np.searchsorted(depths, max_depth, side="left") + 1]
        expected_grid = {}
        for name, lower, upper in (("lat", "south", "north"), ("lon", "west", "east")):
            values = (
                axes[name] if name == "lat" else np.sort((axes[name] + 180) % 360 - 180)
            )
            # Modulo conversion can move an exact endpoint by a few ULPs.
            tolerance = 1e-10
            if (
                settings[lower] < values[0] - tolerance
                or settings[upper] > values[-1] + tolerance
            ):
                raise ValueError(
                    f"Requested {name} bounds exceed the HYCOM source grid."
                )
            lower_value = max(settings[lower], values[0])
            upper_value = min(settings[upper], values[-1])
            left = int(np.searchsorted(values, lower_value, side="right")) - 1
            right = int(np.searchsorted(values, upper_value, side="left"))
            expected_grid[name] = values[left : right + 1]
            if right <= left:
                raise ValueError(
                    "The requested region needs at least two native latitude/longitude nodes."
                )
        options = {
            "expected_times": expected_times,
            "expected_depths": retained,
            "expected_grid": expected_grid,
            "variables": tuple(COMPONENTS.values()),
        }
        if output_path.is_file() and not force:
            try:
                with xr.open_dataset(output_path, engine="netcdf4") as source:
                    if source.attrs.get("request_signature") != signature or any(
                        source.attrs.get(name) != value
                        for name, value in signature_data.items()
                    ):
                        raise ValueError("Cached output belongs to another request.")
                    return _validate(source, **options)
            except (OSError, ValueError, RuntimeError):
                pass
        pieces = []
        left = start
        while left < end:
            right = min(left + pd.Timedelta(days=1), end)
            chunk_times = pd.date_range(left, right, freq="3h")
            chunk_path = cache / f"{left:%Y%m%dT%H%M}_{right:%Y%m%dT%H%M}.nc"
            validation = {**options, "expected_times": chunk_times}
            chunk = None
            if chunk_path.is_file() and not force:
                try:
                    with xr.open_dataset(chunk_path, engine="netcdf4") as source:
                        if source.attrs.get("request_signature") != signature or any(
                            source.attrs.get(name) != value
                            for name, value in signature_data.items()
                        ):
                            raise ValueError("Daily cache belongs to another request.")
                        chunk = _validate(source, **validation)
                except (OSError, ValueError, RuntimeError):
                    pass
            if chunk is None:
                fragments = {name: [] for name in COMPONENTS.values()}
                for url, variables in zip(
                    archive["source_urls"], archive["variable_groups"], strict=True
                ):
                    source_times = chunk_times.intersection(axes["source_times"][url])
                    if source_times.empty:
                        continue

                    maximum_times = int(chunk_hours // 3) + 1
                    for index in range(0, len(source_times), maximum_times):
                        packet_times = source_times[index : index + maximum_times]
                        fragment = _fetch_fragment(
                            session,
                            url,
                            variables=variables,
                            times=packet_times,
                            depths=depths,
                            retained=retained,
                            grid=expected_grid,
                            cache=cache,
                            signature_data=signature_data,
                            signature=signature,
                            force=force,
                            method=method,
                            request_options=request_options,
                            source_axes=axes.get("source_axes", {}).get(url),
                        )
                        for variable in variables:
                            fragments[variable].append(fragment[[variable]])
                components = []
                for variable, records in fragments.items():
                    component = xr.concat(
                        records,
                        dim="time",
                        data_vars="minimal",
                        coords="minimal",
                        compat="equals",
                        join="exact",
                        combine_attrs="drop_conflicts",
                    )
                    # Experiments are ordered oldest first; a newer archive is
                    # authoritative at a shared timestamp, without filling gaps.
                    keep = ~pd.Index(component.time.values).duplicated(keep="last")
                    component = component.isel(time=np.flatnonzero(keep)).sortby("time")
                    components.append(
                        _validate(component, **{**validation, "variables": [variable]})
                    )
                chunk = xr.merge(
                    components,
                    join="exact",
                    compat="equals",
                    combine_attrs="drop_conflicts",
                )
                chunk.attrs.update(
                    signature_data,
                    request_signature=signature,
                )
                _save(chunk, chunk_path)
            pieces.append(chunk)
            left = right
    combined = xr.concat(
        pieces,
        dim="time",
        data_vars="minimal",
        coords="minimal",
        compat="equals",
        join="exact",
        combine_attrs="drop_conflicts",
    )
    keep = ~pd.Index(combined.time.values).duplicated(keep="last")
    combined = combined.isel(time=np.flatnonzero(keep))
    combined = _validate(combined, **options)
    combined.attrs.update(
        signature_data,
        request_signature=signature,
        time_coverage_start=start.isoformat(),
        time_coverage_end=end.isoformat(),
        requested_start=start.isoformat(),
        requested_end=end.isoformat(),
    )
    _save(combined, output_path)
    return combined
