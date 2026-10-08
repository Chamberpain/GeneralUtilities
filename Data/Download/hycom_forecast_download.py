"""Regional HYCOM forecast subsets on native model depths.

No datasets are opened until a public function is called. Forecast runs are
stitched locally, preferring the newer run when timestamps overlap.
"""

import re
import time
import warnings
import xml.etree.ElementTree as ET

import numpy as np
import pandas as pd
import requests
import xarray as xr

from ._subset_download_utils import read_netcdf_bytes, request_bytes


# ============================================================
# TIME FORMAT
# ============================================================

def return_tds_dataset(tracer=False):
    base_html = "https://tds.hycom.org/thredds/dodsC/"
    if tracer:
        url = f"{base_html}FMRC_ESPC-D-V02_ts3z/FMRC_ESPC-D-V02_ts3z_best.ncd"
        print('I am opening the tracer dataset')
    else:
        url = f"{base_html}FMRC_ESPC-D-V02_uv3z/FMRC_ESPC-D-V02_uv3z_best.ncd"
        print('I am opening the uv dataset')
    session = requests.Session()
    session.headers.update({
        "User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36"
    })

    # The live remote Dataset keeps the PyDAP session available until closed.
    for attempt in range(2):
        try:
            return xr.open_dataset(
                url, engine="pydap", drop_variables=["tau"],
                backend_kwargs={"session": session, "timeout": 500},
            )
        except Exception:
            if attempt:
                session.close()
                raise
            print("I encountered an error opening the dataset, trying again")


def NCSS_dataset(
    west, east, south, north, start, end,
    tracer=False, max_depth=500.0, time_stride=1, horiz_stride=1,
    max_attempts=5, timeout=500, retry_wait=2,
) -> xr.Dataset:
    """Download U/V (or T/S with tracer=True) on every native depth <= max_depth.

    Stitch individual forecast runs; keep the newer run at duplicate times.
    Retry failed requests, but never silently skip a failed depth/run segment.
    Warn on partial time coverage. time_stride is applied within each segment.
    timeout is the HTTP read timeout, not a total download time limit.
    """
    try:
        import netCDF4  # noqa: F401 -- check the backend before downloading
    except ImportError as exc:
        raise ImportError("Install netCDF4: conda install -c conda-forge netcdf4") from exc

    start = pd.to_datetime(start, utc=True)
    end = pd.to_datetime(end, utc=True)
    if pd.isna(start) or pd.isna(end) or end <= start:
        raise ValueError("Provide valid dates with end after start.")
    if not np.isfinite([west, east, south, north, max_depth, timeout, retry_wait]).all():
        raise ValueError("Bounds, maximum depth, timeout, and retry_wait must be finite.")
    if not -90 <= south < north <= 90:
        raise ValueError("Require -90 <= south < north <= 90.")
    if east < west:
        east += 360
    if not 0 < east - west <= 360:
        raise ValueError("Longitude bounds must span more than 0 and at most 360 degrees.")
    if max_depth < 0 or timeout <= 0 or retry_wait < 0:
        raise ValueError("Require max_depth >= 0, timeout > 0, and retry_wait >= 0.")
    for name, value in (("time_stride", time_stride), ("horiz_stride", horiz_stride),
                        ("max_attempts", max_attempts)):
        if isinstance(value, bool) or not isinstance(value, (int, np.integer)) or value < 1:
            raise ValueError(f"{name} must be a positive integer.")

    product = "ts3z" if tracer else "uv3z"
    variables = ["water_temp", "salinity"] if tracer else ["water_u", "water_v"]
    prefix = f"FMRC_ESPC-D-V02_{product}"
    catalog_url = f"https://tds.hycom.org/thredds/catalog/{prefix}/runs/catalog.xml"
    ncss_base = "https://ncss.hycom.org/thredds/ncss/grid/"

    # Parsing and model-specific validation share the common retry transport.
    def request(url, params=None, *, label="request", reference=None):
        def decode(content):
            if params is None:
                root = ET.fromstring(content)
                if root.tag.rsplit("}", 1)[-1] != "catalog":
                    raise ValueError("Expected a THREDDS XML catalog.")
                return root

            ds = read_netcdf_bytes(content, prefix="hycom_")

            wanted = params["var"].split(",")
            if any(v not in ds.data_vars or ds[v].size == 0 for v in wanted):
                raise ValueError(f"Missing or empty variables; expected {wanted}.")
            if "time" not in ds or ds.time.size == 0:
                raise ValueError("Missing or empty time coordinate.")

            # Point profiles may use a different name for their depth variable.
            z = "depth" if "depth" in ds else next((
                name for name, var in ds.variables.items()
                if var.attrs.get("axis", "").upper() == "Z"
                or var.attrs.get("standard_name", "") == "depth"
                or var.attrs.get("positive", "").lower() == "down"
            ), None)
            if z is None:
                raise ValueError("No depth coordinate found.")
            if z != "depth":
                ds = ds.rename({z: "depth"})
            ds = ds.set_coords("depth")
            depths = np.asarray(ds.depth.values, dtype=float).ravel()
            if not depths.size or not np.isfinite(depths).all():
                raise ValueError("Invalid depth coordinates.")
            units = ds.depth.attrs.get("units", "m").strip().lower()
            if units not in ("m", "meter", "meters", "metre", "metres"):
                raise ValueError(f"Expected depths in meters, received {units!r}.")
            if ds.depth.attrs.get("positive", "down").lower() != "down":
                raise ValueError("Expected positive-down HYCOM depths.")

            if "vertCoord" in params:
                if depths.size != 1 or not np.isclose(
                    depths[0], params["vertCoord"], rtol=0, atol=1e-4,
                ):
                    raise ValueError(f"Wrong depth: requested {params['vertCoord']}, got {depths}.")
                if "depth" not in ds.dims:
                    ds = ds.expand_dims("depth")
                required = {"time", "depth", "lat", "lon"}
                if any(set(ds[v].dims) != required for v in wanted):
                    raise ValueError("Expected each field to have time/depth/lat/lon dimensions.")
                # Discard per-run scalar coordinates that would conflict when stitching.
                ds = ds[wanted].reset_coords(drop=True)
                times = ds.time.to_index()
                if times.hasnans or not times.is_unique or not times.is_monotonic_increasing:
                    raise ValueError("Time coordinate is invalid, duplicated, or unsorted.")
                if reference is not None:
                    for coord in ("time", "lat", "lon"):
                        if not ds[coord].equals(reference[coord]):
                            raise ValueError(f"{coord} differs between depths for this run segment.")
            return ds


        def report(error, attempt, delay):
            print(f"  {label}: {error}\n  Retry {attempt + 2}/{max_attempts} in {delay:g}s.", flush=True)

        return request_bytes(
            session, url, params=params, timeout=timeout, max_attempts=max_attempts,
            retry_wait=retry_wait, reader=decode, connect_timeout=30, retry_cap=60,
            failure_message=f"{label} failed after {max_attempts} attempts.",
            retry_reporter=report,
        )

    print(f"HYCOM {', '.join(variables)} | {start} -> {end} | 0-{max_depth:g} m")
    with requests.Session() as session:
        session.headers.update({"User-Agent": "HYCOM-NCSS-Python/1.0"})

        # Discover the run catalog once, retaining the actual advertised paths.
        root = request(catalog_url, label="Run catalog")
        run_urls = {}
        for node in root.iter():
            path = node.get("urlPath", "")
            match = re.search(r"_RUN_(\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}Z)", path)
            if path.startswith(f"{prefix}/runs/") and match:
                run_urls[pd.Timestamp(match[1])] = ncss_base + path
        runs = sorted(run_urls)
        if not runs:
            raise RuntimeError(f"No individual forecast runs found in {catalog_url}.")
        print(f"Found {len(runs)} runs: {runs[0]} -> {runs[-1]}")

        # Use the latest run at/before start, followed by each new run before end.
        prior = [run for run in runs if run <= start]
        selected = prior[-1:] + [run for run in runs if start < run < end]
        if not selected:
            raise ValueError("No run begins before the requested end time.")
        segments = [
            (run, max(start, run), min(end, selected[i + 1]) if i + 1 < len(selected) else end)
            for i, run in enumerate(selected)
        ]

        # A small ocean-point profile returns all native depths; no vertCoord is sent.
        profile = request(run_urls[runs[-1]], {
            "var": variables[0], "longitude": -140, "latitude": 30,
            "time": runs[-1].isoformat(), "accept": "netCDF4",
        }, label="Depth discovery")
        depths = np.unique(np.asarray(profile.depth.values, dtype=float))
        depths = depths[(depths >= 0) & (depths <= max_depth)]
        if not depths.size:
            raise ValueError(f"No native depths between 0 and {max_depth:g} m.")
        print(f"Downloading {len(depths)} depths across {len(segments)} run segments:\n{depths}")

        params = {
            "var": ",".join(variables), "west": west, "east": east,
            "south": south, "north": north, "timeStride": int(time_stride),
            "horizStride": int(horiz_stride), "accept": "netCDF4", "addLatLon": "true",
        }
        concat_options = dict(
            data_vars="minimal", coords="minimal", compat="equals",
            join="exact", combine_attrs="drop_conflicts",
        )
        layers, references = [], {}
        for i, depth in enumerate(depths, 1):
            print(f"[{i}/{len(depths)}] Depth {depth:g} m", flush=True)
            pieces = []
            for run, left, right in segments:
                part = request(run_urls[run], {
                    **params, "vertCoord": float(depth),
                    "time_start": left.isoformat(), "time_end": right.isoformat(),
                }, label=f"{depth:g} m / run {run}", reference=references.get(run))
                references.setdefault(run, part[["time", "lat", "lon"]])
                pieces.append(part)

            layer = xr.concat(pieces, dim="time", **concat_options)
            # Deduplicate BEFORE sorting: pieces are ordered from older to newer runs.
            keep = ~pd.Index(layer.time.values).duplicated(keep="last")
            layer = layer.isel(time=np.flatnonzero(keep)).sortby("time")
            layer = layer.sel(time=slice(start.tz_localize(None), end.tz_localize(None)))
            if not layer.sizes["time"]:
                raise ValueError("No returned times fall inside the requested interval.")
            layers.append(layer)
            print(f"  Complete: {layer.sizes['time']} times", flush=True)

    dataset = xr.concat(layers, dim="depth", **concat_options).sortby("depth")
    dataset = dataset.transpose("time", "depth", "lat", "lon")
    if not np.allclose(dataset.depth.values, depths, rtol=0, atol=1e-4):
        raise RuntimeError("The combined depth grid does not match the requested native grid.")
    actual_start, actual_end = pd.to_datetime(dataset.time.values[[0, -1]], utc=True)
    if actual_start > start or actual_end < end:
        warnings.warn(
            f"Partial time coverage. Requested {start} -> {end}; "
            f"received {actual_start} -> {actual_end}. No extrapolation was performed.",
            RuntimeWarning, stacklevel=2,
        )
    dataset.attrs.update(
        download_method="HYCOM NCSS individual runs stitched locally",
        maximum_depth_m=float(max_depth), requested_start=start.isoformat(),
        requested_end=end.isoformat(), time_coverage_start=actual_start.isoformat(),
        time_coverage_end=actual_end.isoformat(),
        source_runs=",".join(run.isoformat() for run in selected),
    )
    print(f"Done: {dict(dataset.sizes)}\nReceived: {actual_start} -> {actual_end}")
    return dataset
