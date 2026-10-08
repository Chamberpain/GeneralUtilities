"""Regional Copernicus Marine current subsets independent of project consumers.

Credentials are passed directly to the optional SDK. Publication occurs only
after data loading, consumer validation, and requested coverage checks succeed.
"""

from __future__ import annotations

import importlib
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np
import pandas as pd
import xarray as xr

DEFAULT_DATASET_ID = "cmems_mod_glo_phy-cur_anfc_0.083deg_PT6H-i"


def validate_copernicus_request(start, end, *, max_depth, dataset_id, output_path=None):
    """Validate inexpensive options before callers compute their model region."""
    start, end = pd.to_datetime(start, utc=True), pd.to_datetime(end, utc=True)
    if pd.isna(start) or pd.isna(end) or end <= start:
        raise ValueError("Provide valid dates with end after start.")
    if (
        isinstance(max_depth, bool)
        or not np.isscalar(max_depth)
        or not np.isfinite(max_depth)
        or max_depth < 0
    ):
        raise ValueError("max_depth must be a finite nonnegative depth in meters.")
    if not isinstance(dataset_id, str) or not dataset_id.strip():
        raise ValueError("dataset_id must be a nonempty Copernicus dataset ID.")
    destination = Path(output_path).expanduser() if output_path is not None else None
    if destination is not None and destination.suffix.lower() != ".nc":
        raise ValueError("output_path must name a .nc NetCDF file.")
    return start, end, destination


def download_copernicus_data(
    start,
    end,
    *,
    bounds,
    output_path=None,
    dataset_id=DEFAULT_DATASET_ID,
    max_depth=500.0,
    username=None,
    password=None,
    validate_dataset=None,
):
    """Return loaded U/V currents covering inclusive UTC dates and geographic bounds.

    Native grid points outside bounds and the requested maximum depth are
    retained for interpolation. An optional ``validate_dataset(dataset)`` hook
    runs before publication, allowing downstream FieldSet validation without a
    dependency on Parcels or project-specific waypoint generators.
    """
    start, end, destination = validate_copernicus_request(
        start, end, max_depth=max_depth, dataset_id=dataset_id, output_path=output_path,
    )
    if set(bounds) != {"west", "east", "south", "north"}:
        raise ValueError("bounds must contain west/east/south/north.")
    bounds = {name: float(value) for name, value in bounds.items()}
    if not np.isfinite(list(bounds.values())).all() or not (
        -180 <= bounds["west"] < bounds["east"] <= 180
        and -90 <= bounds["south"] < bounds["north"] <= 90
    ):
        raise ValueError("Provide finite geographic bounds without a dateline crossing.")
    try:
        sdk = importlib.import_module("copernicusmarine")
    except ImportError as error:
        raise ImportError(
            "Copernicus downloads require the optional SDK: "
            "python -m pip install copernicusmarine. Then run copernicusmarine login."
        ) from error

    if destination is not None:
        destination.parent.mkdir(parents=True, exist_ok=True)
    options = {
        "dataset_id": dataset_id,
        "variables": ["uo", "vo"],
        "minimum_longitude": bounds["west"],
        "maximum_longitude": bounds["east"],
        "minimum_latitude": bounds["south"],
        "maximum_latitude": bounds["north"],
        "minimum_depth": 0.0,
        "maximum_depth": float(max_depth),
        "vertical_axis": "depth",
        "start_datetime": start.isoformat(),
        "end_datetime": end.isoformat(),
        "coordinates_selection_method": "outside",
        "output_filename": "currents.nc",
        "file_format": "netcdf",
        "overwrite": True,
    }
    options.update(
        {
            key: value
            for key, value in {"username": username, "password": password}.items()
            if value is not None
        }
    )
    # Stage beside the destination so replace is atomic on the same filesystem.
    with TemporaryDirectory(
        prefix="copernicus_",
        dir=destination.parent if destination is not None else None,
    ) as folder:
        sdk.subset(**options, output_directory=folder)
        temporary = Path(folder) / "currents.nc"
        if not temporary.is_file():
            raise RuntimeError(
                "The Copernicus SDK did not create the requested NetCDF file."
            )
        with xr.open_dataset(temporary) as source:
            dataset = source.load()
        if not {"uo", "vo"}.issubset(dataset.data_vars):
            raise ValueError("Downloaded currents must include uo and vo.")
        if validate_dataset is not None:
            validate_dataset(dataset)
        times = pd.to_datetime(dataset.time.values, utc=True)
        if times.min() > start or times.max() < end:
            raise ValueError(
                "Downloaded currents do not bracket the requested time range. "
                "Choose dates within the dataset's available analysis/forecast record."
            )
        longitude = dataset["longitude" if "longitude" in dataset.coords else "lon"]
        longitude = (longitude + 180.0) % 360.0 - 180.0
        latitude = dataset["latitude" if "latitude" in dataset.coords else "lat"]
        if (
            float(longitude.min()) > bounds["west"]
            or float(longitude.max()) < bounds["east"]
            or float(latitude.min()) > bounds["south"]
            or float(latitude.max()) < bounds["north"]
        ):
            raise ValueError(
                "Downloaded currents do not cover the padded waypoint region."
            )
        if float(dataset.depth.max()) < max_depth:
            raise ValueError(
                "Downloaded currents do not bracket the requested maximum depth."
            )
        if destination is not None:
            temporary.replace(destination)
    return dataset
