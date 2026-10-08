"""WCOFS archive/forecast downloading with physical-grid normalization."""
from __future__ import annotations
import datetime
from contextlib import nullcontext
import hashlib
from pathlib import Path
import pickle
import numpy as np
import pandas as pd
import xarray as xr
from .current_archive import NETCDF_LOCK, download_record as download_archive


def _utc(value):
    value = pd.Timestamp(value)
    if pd.isna(value):
        raise ValueError("A finite UTC date is required")
    return value.tz_localize(None) if value.tz is None else value.tz_convert("UTC").tz_localize(None)


class WCOFSDownloader:
    dataset_description = "WGOFS"
    hours_list = list(range(0, 73, 3))
    time_step = datetime.timedelta(hours=3)
    url_template = (
        "https://opendap.co-ops.nos.noaa.gov/thredds/dodsC/"
        "NOAA/WCOFS/MODELS/{year:04d}/{month:02d}/{day:02d}/"
        "wcofs.t03z.{date}.regulargrid.{kind}{hour:03d}.nc"
    )
    roms_url_template = (
        "https://www.ncei.noaa.gov/thredds/dodsC/model-wcofs-files/"
        "{year:04d}/{month:02d}/{archive_filename}"
    )

    @classmethod
    def source_url(
        cls, valid_time, *, url_template=None, cycle=None, source="regulargrid"
    ):
        """Map an exact 3-hour valid time to a daily 03Z nowcast.

        n003 is cycle minus 21 hours, n024 is the cycle time. For a forecast,
        explicitly supply its 03Z cycle; forecasts never silently replace a
        historical nowcast. Template fields refer to the cycle date.
        """
        if source not in {"regulargrid", "roms"}:
            raise ValueError("WCOFS source must be regulargrid or roms")
        valid_time = _utc(valid_time)
        if (
            valid_time.minute
            or valid_time.second
            or valid_time.microsecond
            or valid_time.nanosecond
            or valid_time.hour % 3
        ):
            raise ValueError(
                "WCOFS regular-grid valid times must lie on the 3-hour UTC grid"
            )
        if cycle is None:
            cycle = valid_time.normalize() + pd.Timedelta(hours=3)
            if valid_time > cycle:
                cycle += pd.Timedelta(days=1)
            hour = int((valid_time - cycle).total_seconds() / 3600) + 24
            kind = "n"
        else:
            cycle = _utc(cycle)
            if (
                cycle.hour != 3
                or cycle.minute
                or cycle.second
                or cycle.microsecond
                or cycle.nanosecond
            ):
                raise ValueError(
                    "A WCOFS forecast cycle must be an exact 03Z UTC timestamp"
                )
            hour = int((valid_time - cycle).total_seconds() / 3600)
            if hour == 0:
                kind, hour = "n", 24
            elif hour < 3 or hour > 72 or hour % 3:
                raise ValueError(
                    "Forecast times must be 3 to 72 hours after the selected cycle"
                )
            else:
                kind = "f"
        product = "regulargrid" if source == "regulargrid" else "fields"
        date = cycle.strftime("%Y%m%d")
        filename = f"wcofs.t03z.{date}.{product}.{kind}{hour:03d}.nc"
        legacy_filename = f"nos.wcofs.{product}.{kind}{hour:03d}.{date}.t03z.nc"
        archive_filename = (
            legacy_filename if cycle < pd.Timestamp("2024-09-10") else filename
        )
        template = url_template or (
            cls.roms_url_template if source == "roms" else cls.url_template
        )
        return template.format(
            year=cycle.year,
            month=cycle.month,
            day=cycle.day,
            date=cycle.strftime("%Y%m%d"),
            cycle=cycle.strftime("%Y-%m-%dT%H:%M:%S"),
            kind=kind,
            hour=hour,
            valid_time=valid_time.strftime("%Y-%m-%dT%H:%M:%S"),
            filename=filename,
            legacy_filename=legacy_filename,
            archive_filename=archive_filename,
        )

    @classmethod
    def get_dataset(
        cls,
        datetime_instance,
        *,
        url_template=None,
        cycle=None,
        engine="pydap",
        source="regulargrid",
    ):
        """Open the requested valid time lazily; the caller must close it."""
        path = cls.source_url(
            datetime_instance, url_template=url_template, cycle=cycle, source=source
        )
        kwargs = {} if engine is None else {"engine": engine}
        try:
            return xr.open_dataset(path, **kwargs)
        except (OSError, ValueError, ImportError) as exc:
            availability = (
                "NOAA's regular-grid archive retains approximately three days; "
                "older dates need a retained regular-grid archive. "
                "Raw ROMS fields require source='roms' and physical preprocessing."
                if source == "regulargrid"
                else "Check the historical raw-field archive URL, naming convention "
                "and availability; no replacement date will be used."
            )
            raise OSError(
                f"Cannot open WCOFS {source} valid time {_utc(datetime_instance)} at {path}. "
                f"{availability}"
            ) from exc

    @staticmethod
    def _normalise(dataset, *, bounds, max_depth):
        """Subset before loading velocity, verifying the source grid and units."""
        dataset = xr.decode_cf(dataset)
        names = {}
        for target, aliases in {
            "U": ("u_eastward", "U"),
            "V": ("v_northward", "V"),
            "time": ("time",),
            "depth": ("Depth", "depth"),
            "lat": ("Latitude", "latitude", "lat"),
            "lon": ("Longitude", "longitude", "lon"),
        }.items():
            matches = [name for name in aliases if name in dataset]
            if not matches:
                if target == "depth":
                    continue
                raise ValueError(
                    f"Missing WCOFS regular-grid {target}; raw ROMS fields need preprocessing"
                )
            names[target] = matches[0]
        lat, lon = dataset[names["lat"]], dataset[names["lon"]]
        if lat.ndim == lon.ndim == 2:
            if lat.dims != lon.dims:
                raise ValueError(
                    "Latitude and longitude must use the same horizontal grid"
                )
            lat_values, lon_values = lat.values, lon.values
            if not np.allclose(
                lat_values, lat_values[:, :1], rtol=0, atol=1e-6
            ) or not np.allclose(lon_values, lon_values[:1, :], rtol=0, atol=1e-6):
                raise ValueError(
                    "Curvilinear coordinates require horizontal regridding"
                )
            lat_dim, lon_dim = lat.dims
            lat_values, lon_values = lat_values[:, 0], lon_values[0, :]
        elif lat.ndim == lon.ndim == 1:
            lat_dim, lon_dim = lat.dims[0], lon.dims[0]
            lat_values, lon_values = lat.values, lon.values
        else:
            raise ValueError(
                "WCOFS latitude/longitude must be one-dimensional or a separable regular grid"
            )
        lon_values = (np.asarray(lon_values, dtype=float) + 180) % 360 - 180
        lat_values = np.asarray(lat_values, dtype=float)
        time = dataset[names["time"]]
        if time.ndim != 1 or not np.issubdtype(time.dtype, np.datetime64):
            raise ValueError(
                "WCOFS time must decode to a one-dimensional UTC datetime axis"
            )
        dim_map = {time.dims[0]: "time", lat_dim: "lat", lon_dim: "lon"}
        if "depth" in names:
            depth = dataset[names["depth"]]
            if depth.ndim != 1 or str(
                depth.attrs.get("units", "")
            ).lower().strip() not in {"m", "meter", "meters", "metre", "metres"}:
                raise ValueError(
                    "WCOFS depths must be a one-dimensional physical metre axis; sigma coordinates need regridding"
                )
            depth_values = np.asarray(depth.values, dtype=float)
            positive = str(depth.attrs.get("positive", "down")).lower()
            if positive == "up":
                depth_values = -depth_values
            elif positive != "down":
                raise ValueError("WCOFS depth positive direction must be up or down")
            if not np.isfinite(depth_values).all() or (depth_values < 0).any():
                raise ValueError(
                    "WCOFS physical depths must be finite and positive down"
                )
            dim_map[depth.dims[0]] = "depth"
        else:
            depth_values = np.asarray([0.0])
        result = xr.Dataset(
            coords={
                "time": ("time", time.values),
                "depth": ("depth", depth_values),
                "lat": ("lat", lat_values),
                "lon": ("lon", lon_values),
            }
        )
        for target in ("U", "V"):
            source = dataset[names[target]]
            unexpected = set(source.dims) - set(dim_map)
            if unexpected:
                raise ValueError(
                    f"{target} has unsupported/staggered dimensions: {sorted(unexpected)}"
                )
            if "depth" in names and dataset[names["depth"]].dims[0] not in source.dims:
                raise ValueError("Current variables must use the physical depth axis")
            # Renaming dimensions does not load remote velocity data.
            dims = tuple(dim_map[name] for name in source.dims)
            value = xr.DataArray(source.variable).rename(
                {name: dim_map[name] for name in source.dims}
            )
            if "depth" not in dims:
                value = value.expand_dims(depth=1)
            if set(value.dims) != {"time", "depth", "lat", "lon"}:
                raise ValueError(
                    f"{target} must span time, depth, latitude and longitude"
                )
            units = str(source.attrs.get("units", "")).lower().strip()
            if units not in {
                "m/s",
                "m s-1",
                "m s^-1",
                "m s**-1",
                "meter second-1",
                "meters second-1",
                "meters/second",
                "metres/second",
            }:
                raise ValueError(
                    f"WCOFS {target} must declare metre-per-second units, got {units!r}"
                )
            result[target] = value.transpose("time", "depth", "lat", "lon")
        for axis in ("depth", "lat", "lon"):
            values = result[axis].values
            if not np.isfinite(values).all() or len(np.unique(values)) != len(values):
                raise ValueError(f"WCOFS {axis} must contain unique finite coordinates")
            result = result.sortby(axis)
        west, east, south, north = map(float, bounds)
        if not all(np.isfinite(bounds)) or not (west < east and south < north):
            raise ValueError("Bounds must be finite west < east and south < north")
        indexers = {}
        for axis, low, high in (("lat", south, north), ("lon", west, east)):
            values = result[axis].values
            if low < values[0] or high > values[-1]:
                raise ValueError(
                    f"Requested {axis} bounds extend beyond WCOFS coverage"
                )
            lower = max(int(np.searchsorted(values, low, side="right")) - 1, 0)
            upper = min(
                int(np.searchsorted(values, high, side="left")) + 1, len(values)
            )
            indexers[axis] = slice(lower, upper)
        if max_depth is not None:
            maximum = abs(float(max_depth))
            if not np.isfinite(maximum):
                raise ValueError("max_depth must be finite")
            selected = np.flatnonzero(result.depth.values <= maximum)
            if not len(selected):
                raise ValueError("No WCOFS depths fall within max_depth")
            indexers["depth"] = selected
        result = result.isel(indexers).load()
        for target, standard_name in (
            ("U", "eastward_sea_water_velocity"),
            ("V", "northward_sea_water_velocity"),
        ):
            result[target].attrs = {"units": "m s-1", "standard_name": standard_name}
        result.depth.attrs = {"units": "m", "positive": "down", "axis": "Z"}
        result.lat.attrs = {
            "units": "degrees_north",
            "axis": "Y",
            "standard_name": "latitude",
        }
        result.lon.attrs = {
            "units": "degrees_east",
            "axis": "X",
            "standard_name": "longitude",
        }
        result.time.attrs = {"standard_name": "time", "axis": "T"}
        return result

    @classmethod
    def download_record(
        cls,
        start,
        end,
        output_dir,
        *,
        workers=4,
        chunk_days=1,
        retries=3,
        worker_index=None,
        max_depth=None,
        bounds=None,
        url_template=None,
        cycle=None,
        engine="pydap",
        source="roms",
        depth_levels=None,
        grid_spacing=0.05,
        preprocess=None,
    ):
        """Download nowcasts in parallel, resuming validated canonical chunks.

        Supply a fixed ``cycle`` only for a particular forecast (up to 72h).
        ``url_template`` also accepts local filesystem templates for retained
        NOAA regular-grid files. Never replace missing times with today's data.
        The default ``source='roms'`` uses NCEI historical nowcasts and performs
        destaggering, rotation and physical depth/geographic interpolation.
        Forecasts require ``source='regulargrid'`` or a retained forecast URL.
        Templates may use ``archive_filename`` (date-correct NOAA naming),
        ``filename`` or ``legacy_filename``. A custom preprocessing callback
        receives dataset, bounds, max_depth, depth_levels and grid_spacing.
        """
        if source not in {"regulargrid", "roms"}:
            raise ValueError("WCOFS source must be regulargrid or roms")
        if source == "regulargrid" and (
            preprocess is not None or depth_levels is not None
        ):
            raise ValueError("ROMS preprocessing options require source='roms'")
        processor_hash = None
        processor_id = None
        if source == "roms" and preprocess is None:
            from .wcofs_roms import regrid_roms

            preprocess = regrid_roms
            processor_id = "HyperNav.Utilities.Data.wcofs_roms:regrid_roms"
            processor_hash = hashlib.sha256(
                Path(__file__).with_name("wcofs_roms.py").read_bytes()
            ).hexdigest()
        if (
            getattr(preprocess, "__module__", "") == "GeneralUtilities.Data.Download.wcofs_roms"
            and getattr(preprocess, "__qualname__", "") == "regrid_roms"
        ):
            processor_id = "HyperNav.Utilities.Data.wcofs_roms:regrid_roms"
        if preprocess is not None and not callable(preprocess):
            raise TypeError("preprocess must be callable")
        bounds = tuple(
            (cls.lllon, cls.urlon, cls.lllat, cls.urlat) if bounds is None else bounds
        )
        if (
            len(bounds) != 4
            or not np.isfinite(bounds).all()
            or not (bounds[0] < bounds[1] and bounds[2] < bounds[3])
        ):
            raise ValueError("Bounds must be finite west < east and south < north")
        max_depth = cls.max_depth if max_depth is None else max_depth
        if isinstance(max_depth, (bool, np.bool_)) or not np.isfinite(float(max_depth)):
            raise ValueError("max_depth must be a finite physical depth")
        max_depth = abs(float(max_depth))
        template = url_template or (
            cls.roms_url_template if source == "roms" else cls.url_template
        )
        cycle = None if cycle is None else _utc(cycle)

        def fetch_chunk(chunk_start, chunk_end):
            times = pd.date_range(
                _utc(chunk_start).ceil("3h"),
                _utc(chunk_end),
                freq="3h",
                inclusive="left",
            )
            if not len(times):
                raise ValueError("Chunk contains no WCOFS 3-hour valid times")
            pieces, paths = [], []
            for valid_time in times:
                path = cls.source_url(
                    valid_time, url_template=template, cycle=cycle, source=source
                )
                # Native NetCDF libraries may share process-level C state.
                with nullcontext() if engine == "pydap" else NETCDF_LOCK:
                    with cls.get_dataset(
                        valid_time,
                        url_template=template,
                        cycle=cycle,
                        engine=engine,
                        source=source,
                    ) as dataset:
                        if source == "roms":
                            part = preprocess(
                                dataset,
                                bounds=bounds,
                                max_depth=abs(float(max_depth)),
                                depth_levels=depth_levels,
                                grid_spacing=grid_spacing,
                            )
                        else:
                            part = cls._normalise(
                                dataset, bounds=bounds, max_depth=max_depth
                            )
                        # Custom processors may retain the open source backend.
                        part.load()
                if (
                    part.sizes["time"] != 1
                    or part.time.values[0] != valid_time.to_datetime64()
                ):
                    raise ValueError(
                        f"WCOFS file {path} contains a different valid time than {valid_time}"
                    )
                pieces.append(part)
                paths.append(path)
            result = xr.concat(pieces, dim="time", join="exact", compat="equals")
            result.attrs.update(
                {
                    "model": "WCOFS",
                    "source_urls": "\n".join(paths),
                    "product": "nowcast" if cycle is None else "selected_cycle",
                }
            )
            return result

        return download_archive(
            fetch_chunk,
            start=start,
            end=end,
            output_dir=output_dir,
            model="WCOFS",
            workers=workers,
            chunk_days=chunk_days,
            retries=retries,
            worker_index=worker_index,
            settings={
                "bounds": list(bounds),
                "max_depth": max_depth,
                "url_template": template,
                "cycle": None if cycle is None else cycle.isoformat(),
                "region": cls.location,
                "source": source,
                "grid_spacing": grid_spacing if source == "roms" else None,
                "depth_levels": None
                if depth_levels is None
                else np.asarray(depth_levels).tolist(),
                "preprocess": None
                if preprocess is None
                else processor_id or f"{getattr(preprocess, '__module__', type(preprocess).__module__)}:"
                f"{getattr(preprocess, '__qualname__', type(preprocess).__qualname__)}",
                "preprocessor_sha256": processor_hash,
            },
        )

    @classmethod
    def download_and_save(cls, start=None, end=None, output_dir=None, **kwargs):
        """Preserve legacy pickle downloads; explicit dates select NetCDF records."""
        if start is not None or end is not None or output_dir is not None:
            if start is None or end is None or output_dir is None:
                raise ValueError(
                    "NetCDF record downloads need start, end and output_dir"
                )
            return cls.download_record(start, end, output_dir, **kwargs)
        if kwargs:
            raise TypeError(
                "Legacy downloads accept no additional options; use download_record"
            )
        now = _utc(pd.Timestamp.now(tz="UTC"))
        cycle = now.normalize() + pd.Timedelta(hours=3)
        if cycle > now:
            cycle -= pd.Timedelta(days=1)
        cls.initialize(cycle, force=True)
        for time in cls.dataset_time:
            filename = Path(cls.file_handler.tmp_file(cls.make_k_filename(time)))
            if filename.exists():
                continue
            with cls.get_dataset(time, cycle=cycle) as source:
                data = cls._normalise(
                    source,
                    bounds=(cls.lllon, cls.urlon, cls.lllat, cls.urlat),
                    max_depth=cls.max_depth,
                )
            filename.parent.mkdir(parents=True, exist_ok=True)
            temporary = filename.with_name(filename.name + ".tmp")
            with temporary.open("wb") as handle:
                pickle.dump(
                    {"u": data.U.values[0], "v": data.V.values[0], "time": time}, handle
                )
            temporary.replace(filename)
