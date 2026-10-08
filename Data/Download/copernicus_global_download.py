"""Copernicus SDK subsets and canonical, resumable current archives.

Subclass with geographic bounds and location, or call get_dataset directly.
SDK authentication remains explicit and managed by copernicusmarine.
"""
import datetime
import importlib
import os
import pickle
import re
import time
from pathlib import Path
import numpy as np
import pandas as pd


class CopernicusDownloader:
    ID = "cmems_mod_glo_phy-cur_anfc_0.083deg_PT6H-i"
    max_depth = 800
    time_step = datetime.timedelta(hours=6)
    hours_list = list(range(0, 25, 6))

    @staticmethod
    def _archive_downloader():
        from .current_archive import download_record
        return download_record

    @classmethod
    def get_dataset(
        cls,
        urlat,
        lllat,
        urlon,
        lllon,
        max_depth,
        ID,
        start_date="2022-06-01",
        end_date=None,
        *,
        variables=None,
        dataset_version=None,
        dataset_part=None,
        credentials_file=None,
        service=None,
    ):
        """Open a lazy SDK subset using saved credentials or SDK environment settings.

        Dates follow the SDK's inclusive selection convention. ``download_record``
        below converts its exclusive stop into an inclusive microsecond timestamp.
        Negative legacy depth limits are accepted as magnitudes in meters.
        """
        depth = abs(float(max_depth))
        if not np.isfinite(depth):
            raise ValueError("max_depth must be a finite depth in meters")
        try:
            sdk = importlib.import_module("copernicusmarine")
        except ImportError as error:
            raise ImportError(
                "Copernicus downloads require copernicusmarine"
            ) from error
        if end_date is None:
            end_date = (datetime.date.today() + datetime.timedelta(days=15)).isoformat()
        options = {
            "dataset_id": ID,
            "minimum_longitude": lllon,
            "maximum_longitude": urlon,
            "minimum_latitude": lllat,
            "maximum_latitude": urlat,
            "start_datetime": start_date,
            "end_datetime": end_date,
            "variables": list(variables or ("uo", "vo")),
            "minimum_depth": 0,
            "maximum_depth": depth,
            "vertical_axis": "depth",
            "coordinates_selection_method": "outside",
        }
        options.update(
            {
                key: value
                for key, value in {
                    "dataset_version": dataset_version,
                    "dataset_part": dataset_part,
                    "credentials_file": credentials_file,
                    "service": service,
                }.items()
                if value is not None
            }
        )
        return sdk.open_dataset(**options)

    @classmethod
    def download_record(
        cls,
        start,
        end,
        output_dir,
        workers=4,
        chunk_days=10,
        retries=3,
        worker_index=None,
        *,
        dataset_id=None,
        max_depth=None,
        time_step=None,
        variables=None,
        dataset_version=None,
        dataset_part=None,
        credentials_file=None,
        service=None,
    ):
        """Save resumable chunks of the half-open UTC interval ``[start, end)``.

        Chunks contain SI ``U,V`` (and ``W`` only if requested) on
        ``time,depth,lat,lon`` with physical positive-down depth in meters.
        ``variables`` maps canonical components to source names, for example
        ``{"U": "uo", "V": "vo", "W": "wo"}``. The default 6-hourly
        operational product starts in June 2022; use a reanalysis dataset ID
        such as ``cmems_mod_glo_phy_my_0.083deg_P1D-m`` for earlier records.
        The SDK's inclusive stop is one microsecond before each exclusive stop.
        Fixed cadence is read from the dataset ID (e.g. PT6H or P1D). Custom
        IDs must supply ``time_step`` as a timedelta, duration string or seconds.
        The actual SDK time coordinates define its phase (daily data may be
        referenced to midnight or noon). Missing native timestamps fail a chunk.
        """
        archive_download_record = cls._archive_downloader()

        components = dict({"U": "uo", "V": "vo"} if variables is None else variables)
        if not {"U", "V"}.issubset(components) or set(components) - {"U", "V", "W"}:
            raise ValueError("variables must map U,V and optionally W to source names")
        if any(not isinstance(name, str) or not name for name in components.values()):
            raise ValueError("Source variable names must be nonempty strings")
        if len(set(components.values())) != len(components):
            raise ValueError("Source velocity variable names must be distinct")
        selected_id = dataset_id or cls.ID
        if time_step is None:
            match = re.search(r"_P(T)?(\d+)([HD])(?:-|_|$)", selected_id)
            if match is None or (bool(match[1]), match[3]) not in (
                (True, "H"),
                (False, "D"),
            ):
                raise ValueError(
                    "Custom dataset IDs require an explicit fixed time_step"
                )
            cadence = (
                pd.Timedelta(hours=int(match[2]))
                if match[3] == "H"
                else pd.Timedelta(days=int(match[2]))
            )
        else:
            if isinstance(time_step, bool):
                raise ValueError("time_step must be a positive fixed duration")
            cadence = (
                pd.Timedelta(seconds=time_step)
                if isinstance(time_step, (int, float))
                else pd.Timedelta(time_step)
            )
        if pd.isna(cadence) or cadence <= pd.Timedelta(0):
            raise ValueError("time_step must be a positive fixed duration")
        depth = abs(float(cls.max_depth if max_depth is None else max_depth))
        if not np.isfinite(depth):
            raise ValueError("max_depth must be finite")
        request = {
            "dataset_id": selected_id,
            "dataset_version": dataset_version,
            "dataset_part": dataset_part,
            "service": service,
            "bounds": [cls.lllon, cls.urlon, cls.lllat, cls.urlat],
            "max_depth": depth,
            "variables": components,
            "time_step_seconds": cadence.total_seconds(),
        }

        def fetch_chunk(chunk_start, chunk_end):
            # Python datetimes/SDK requests have microsecond precision. Clip again
            # locally so even SDK/server boundary rounding cannot duplicate times.
            stop = chunk_end - datetime.timedelta(microseconds=1)
            source = cls.get_dataset(
                cls.urlat,
                cls.lllat,
                cls.urlon,
                cls.lllon,
                depth,
                selected_id,
                chunk_start.isoformat(),
                stop.isoformat(),
                variables=components.values(),
                dataset_version=dataset_version,
                dataset_part=dataset_part,
                credentials_file=credentials_file,
                service=service,
            )
            try:
                selected = source[list(components.values())].rename(
                    {value: key for key, value in components.items()}
                    | {"latitude": "lat", "longitude": "lon"}
                )
                values = selected.time.values
                if len(values) == 0:
                    raise ValueError("Copernicus returned no time coordinates")
                native_phase = (
                    int(np.asarray(values).astype("datetime64[ns]").astype("int64")[0])
                    % cadence.value
                )
                start_value = np.datetime64(chunk_start.value, "ns")
                end_value = np.datetime64(chunk_end.value, "ns")
                selected = selected.isel(
                    time=(values >= start_value) & (values < end_value)
                )
                selected = selected.sortby("time")
                observed = selected.time.values.astype("datetime64[ns]").astype("int64")
                interval_start, interval_end = chunk_start.value, chunk_end.value
                first_expected = (
                    -(-(interval_start - native_phase) // cadence.value) * cadence.value
                    + native_phase
                )
                expected_count = max(
                    0, (interval_end - 1 - first_expected) // cadence.value + 1
                )
                if (
                    len(observed) != expected_count
                    or not len(observed)
                    or observed[0] != first_expected
                    or np.any(np.diff(observed) != cadence.value)
                ):
                    raise ValueError(
                        f"Incomplete Copernicus time coverage: expected {expected_count} "
                        f"records at {cadence}, received {len(observed)}"
                    )
                selected = selected.transpose("time", "depth", "lat", "lon").load()
                scales = {
                    "m s-1": 1.0,
                    "m/s": 1.0,
                    "m s^-1": 1.0,
                    "meter second-1": 1.0,
                    "meters second-1": 1.0,
                    "meters/second": 1.0,
                    "cm s-1": 0.01,
                    "cm/s": 0.01,
                    "mm s-1": 0.001,
                    "mm/s": 0.001,
                }
                for component in components:
                    attrs = dict(selected[component].attrs)
                    units = attrs.get("units", "").strip().lower()
                    if units not in scales:
                        raise ValueError(
                            f"{component} has unsupported velocity units {units!r}"
                        )
                    factor = scales[units]
                    if component == "W":
                        direction = attrs.get("standard_name", "")
                        if direction == "upward_sea_water_velocity":
                            factor *= -1
                        elif direction != "downward_sea_water_velocity":
                            raise ValueError(
                                "W must declare upward/downward_sea_water_velocity"
                            )
                    selected[component] = selected[component] * factor
                    selected[component].attrs = attrs | {
                        "units": "m s-1",
                        "standard_name": {
                            "U": "eastward_sea_water_velocity",
                            "V": "northward_sea_water_velocity",
                            "W": "downward_sea_water_velocity",
                        }[component],
                    }
                if "W" in selected:
                    selected.W.attrs["positive"] = "down"
                depth_attrs = dict(selected.depth.attrs)
                if depth_attrs.get("units", "").lower() not in ("m", "meter", "meters"):
                    raise ValueError("Depth must be expressed in meters")
                if depth_attrs.get("positive", "down") != "down":
                    raise ValueError("SDK must return a positive-down depth axis")
                selected.depth.attrs = depth_attrs | {
                    "standard_name": "depth",
                    "units": "m",
                    "positive": "down",
                    "axis": "Z",
                }
                selected.lat.attrs = dict(selected.lat.attrs) | {
                    "standard_name": "latitude",
                    "units": "degrees_north",
                    "axis": "Y",
                }
                selected = selected.assign_coords(
                    lon=(
                        selected.lon.dims,
                        (np.asarray(selected.lon.values, dtype=float) + 180) % 360
                        - 180,
                        dict(selected.lon.attrs),
                    )
                ).sortby(["depth", "lat", "lon"])
                selected.lon.attrs = dict(selected.lon.attrs) | {
                    "standard_name": "longitude",
                    "units": "degrees_east",
                    "axis": "X",
                }
                selected.time.attrs = dict(selected.time.attrs) | {
                    "standard_name": "time",
                    "axis": "T",
                }
                selected.attrs = dict(selected.attrs) | {
                    "native_time_step_seconds": cadence.total_seconds(),
                    "native_time_phase_seconds": native_phase / 1e9,
                }
                return selected
            finally:
                source.close()

        return archive_download_record(
            fetch_chunk,
            start=start,
            end=end,
            output_dir=output_dir,
            model="Copernicus",
            settings=request,
            workers=workers,
            chunk_days=chunk_days,
            retries=retries,
            worker_index=worker_index,
        )

    @classmethod
    def _download_legacy(cls, start_index=0, *, refresh=False, retries=3):
        cls.initialize()
        indices = cls.dataset_time.return_time_list()
        for index in range(max(0, start_index), len(indices) - 1):
            filename = cls.make_k_filename(index)
            if os.path.isfile(filename) and not refresh:
                continue
            for attempt in range(retries):
                try:
                    dataset = cls.get_dataset(
                        cls.urlat,
                        cls.lllat,
                        cls.urlon,
                        cls.lllon,
                        cls.max_depth,
                        cls.ID,
                        cls.dataset_time[indices[index]].isoformat(),
                        (
                            cls.dataset_time[indices[index + 1]] - cls.time_step
                        ).isoformat(),
                    )
                    try:
                        payload = {
                            "u": dataset.uo.values,
                            "v": dataset.vo.values,
                            "time": dataset.time.values.astype(
                                "datetime64[ns]"
                            ).tolist(),
                        }
                    finally:
                        dataset.close()
                    os.makedirs(os.path.dirname(filename), exist_ok=True)
                    with open(filename, "wb") as handle:
                        pickle.dump(payload, handle)
                    break
                except Exception:
                    if attempt + 1 == retries:
                        raise
                    time.sleep(min(2**attempt, 8))

    @classmethod
    def download_and_save(cls):
        """Retain the original Base pickle format for existing callers."""
        cls._download_legacy()

    @classmethod
    def download_recent(cls):
        cls.initialize()
        cls._download_legacy(
            len(cls.dataset_time.return_time_list()) - 10, refresh=True
        )
