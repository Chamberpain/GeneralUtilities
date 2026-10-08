"""Shared HYCOM/ESPC metadata, regional subsets and resumable pickle downloads.

Subclass HYCOMDownloader with regional bounds, max_depth and a file_handler.
Importing this module neither opens a remote dataset nor initializes a cache.
"""
from __future__ import annotations
import datetime
import json
import os
import pickle
import re
import warnings
from socket import timeout as SocketTimeout
from typing import Any, Dict, Iterable, Mapping, Optional, Sequence, Tuple
from urllib.error import HTTPError
from urllib.parse import urlsplit
import numpy as np
import pandas as pd
import xarray as xr
from requests.exceptions import RequestException
from GeneralUtilities.Compute.list import DepthList, LatList, LonList, TimeList

class CurrentTimeList(TimeList):
    """Legacy current-record chunk/time indexing without simulation dependencies."""
    def return_time_list(self):
        return list(range(len(self))[::100]) + [len(self)-1]

    def get_file_indexes(self, start_time, end_time):
        assert start_time < end_time
        assert start_time > min(self)
        assert end_time < max(self)
        keys = TimeList([self[x] for x in self.return_time_list()])
        start = keys.find_nearest(start_time, idx=True)
        end = keys.find_nearest(end_time, idx=True)
        return np.arange(start-1, end+2).tolist()

# PyDAP is optional until a source is actually opened.
def open_url(*args, **kwargs):
    from pydap.client import open_url as open_source
    return open_source(*args, **kwargs)

try:
    from pydap.exceptions import DapError
except ImportError:
    DapError = RuntimeError


class HYCOMDownloader:
    time_list_class = CurrentTimeList

    @classmethod
    def make_k_filename(cls, k):
        return cls.file_handler.tmp_file(cls.dataset_description + "_" + cls.location + "_data/" + str(k))

    dataset_description = "ESPC_D_V02"
    model_name = "ESPC-D-V02"
    hours_list = np.arange(0, 25, 3).tolist()
    time_step = datetime.timedelta(hours=3)

    # Current operational HYCOM.org FMRC aggregations.  The URL without the
    # trailing ``.html`` is the OPeNDAP endpoint consumed by pydap.
    base_html = "https://tds.hycom.org/thredds/dodsC/"
    velocity_url = base_html + "FMRC_ESPC-D-V02_uv3z/" + "FMRC_ESPC-D-V02_uv3z_best.ncd"
    tracer_url = base_html + "FMRC_ESPC-D-V02_ts3z/" + "FMRC_ESPC-D-V02_ts3z_best.ncd"

    # Kept as a short identifier for callers that inspect/use ``ID``.
    ID = "ESPC-D-V02_uv3z_best"

    # Old IDs are mapped to the new global source so downstream code that
    # still passes one explicitly does not immediately fail.
    _source_aliases = {
        "ESPC-D-V02_uv3z_best": velocity_url,
        "ESPC-D-V02_ts3z_best": tracer_url,
        "HYCOM_reg1_latest3d": velocity_url,
        "HYCOM_reg6_latest3d": velocity_url,
        "HYCOM_reg7_latest3d": velocity_url,
        "HYCOM_reg17_latest3d": velocity_url,
    }
    _legacy_ids = {
        "HYCOM_reg1_latest3d",
        "HYCOM_reg6_latest3d",
        "HYCOM_reg7_latest3d",
        "HYCOM_reg17_latest3d",
    }

    _axis_aliases = {
        "time": ("time",),
        "lat": ("lat", "latitude"),
        "lon": ("lon", "longitude"),
        "depth": ("depth",),
    }

    # Shared across every regional subclass to avoid opening the same global
    # endpoint and downloading the same coordinate axes six times at import.
    _dataset_cache: Dict[str, Any] = {}
    _metadata_cache: Dict[str, Dict[str, Any]] = {}

    dataset = None
    dataset_time = None
    lats = None
    lons = None
    depths = None
    lllon_idx = None
    urlon_idx = None
    lllat_idx = None
    urlat_idx = None
    units = "m/s"
    ref_date = None
    _initialized = False

    @classmethod
    def _resolve_source(cls, source: Optional[str]) -> str:
        source = source or cls.ID
        if source in cls._legacy_ids:
            warnings.warn(
                "%s is retired; using the current ESPC-D-V02 global "
                "velocity feed instead." % source,
                DeprecationWarning,
                stacklevel=3,
            )
        if source in cls._source_aliases:
            return cls._source_aliases[source]
        if source.startswith("http://") or source.startswith("https://"):
            return source
        # Preserve the old ``base_html + ID`` behavior for custom paths.
        return cls.base_html + source.lstrip("/")

    @staticmethod
    def _open_url_compat(url: str, request_timeout: int = 90):
        """Open an OPeNDAP URL across old and new pydap call signatures."""

        try:
            return open_url(
                url,
                timeout=request_timeout,
                output_grid=True,
                protocol="dap2",
            )
        except TypeError:
            try:
                return open_url(url, timeout=request_timeout, output_grid=True)
            except TypeError:
                try:
                    return open_url(url, timeout=request_timeout)
                except TypeError:
                    return open_url(url)

    @classmethod
    def get_dataset(
        cls,
        ID: Optional[str] = None,
        force: bool = False,
        request_timeout: int = 90,
    ):
        """Return a cached pydap dataset for an alias, path, or full URL."""

        url = cls._resolve_source(ID)
        if force:
            HYCOMDownloader._dataset_cache.pop(url, None)
            HYCOMDownloader._metadata_cache.pop(url, None)
        if url not in HYCOMDownloader._dataset_cache:
            HYCOMDownloader._dataset_cache[url] = cls._open_url_compat(
                url, request_timeout=request_timeout
            )
        return HYCOMDownloader._dataset_cache[url]

    @classmethod
    def get_tracer_dataset(cls, force: bool = False):
        return cls.get_dataset(cls.tracer_url, force=force)

    @classmethod
    def archive_url(cls, variable_group: str, year: int) -> str:
        """Build a current ESPC-D-V02 archive URL for advanced callers.

        ``variable_group`` must be one of ``u3z``, ``v3z``, ``t3z`` or
        ``s3z``.  The rolling FMRC feeds above are used by the normal class
        methods; the yearly archive endpoints are exposed here for callers
        that need records older than the rolling best-time-series window.
        """

        valid = {"u3z", "v3z", "t3z", "s3z"}
        if variable_group not in valid:
            raise ValueError(
                "variable_group must be one of %s" % ", ".join(sorted(valid))
            )
        if year < 2024:
            raise ValueError("ESPC-D-V02 archives begin in 2024")
        return "%sESPC-D-V02/%s/%d" % (cls.base_html, variable_group, year)

    @staticmethod
    def _dataset_keys(dataset: Any) -> set:
        try:
            return set(dataset.keys())
        except Exception:
            return set()

    @classmethod
    def _find_variable(cls, dataset: Any, aliases: Sequence[str]) -> str:
        keys = cls._dataset_keys(dataset)
        for name in aliases:
            if not keys or name in keys:
                try:
                    dataset[name]
                    return name
                except (KeyError, TypeError):
                    continue
        raise KeyError(
            "None of the expected variables %r exist in the HYCOM dataset"
            % (tuple(aliases),)
        )

    @staticmethod
    def _attributes(variable: Any) -> Mapping[str, Any]:
        attrs = getattr(variable, "attributes", {})
        return attrs if isinstance(attrs, Mapping) else {}

    @staticmethod
    def _unwrap_data(value: Any, variable_name: Optional[str] = None) -> np.ndarray:
        """Convert pydap BaseType/GridType results to a NumPy array."""

        candidate = value
        if variable_name is not None:
            try:
                candidate = value[variable_name]
            except (KeyError, TypeError, IndexError):
                candidate = value
        if not isinstance(candidate, (np.ndarray, np.ma.MaskedArray)) and hasattr(
            candidate, "data"
        ):
            candidate = candidate.data
        return np.asanyarray(candidate)

    @classmethod
    def _read_axis(cls, dataset: Any, name: str) -> Tuple[str, np.ndarray]:
        variable_name = cls._find_variable(dataset, cls._axis_aliases[name])
        result = dataset[variable_name][:]
        array = cls._unwrap_data(result, variable_name).squeeze()
        return variable_name, np.asarray(array, dtype=float)

    @staticmethod
    def _parse_datetime(text: str) -> datetime.datetime:
        cleaned = text.strip()
        cleaned = re.sub(r"\s+(UTC|GMT)$", "", cleaned, flags=re.IGNORECASE)
        if cleaned.endswith("Z"):
            cleaned = cleaned[:-1] + "+00:00"

        parsed: Optional[datetime.datetime] = None
        try:
            parsed = datetime.datetime.fromisoformat(cleaned)
        except (ValueError, AttributeError):
            pass

        if parsed is None:
            formats = (
                "%Y-%m-%d %H:%M:%S.%f",
                "%Y-%m-%d %H:%M:%S",
                "%Y-%m-%dT%H:%M:%S.%f",
                "%Y-%m-%dT%H:%M:%S",
                "%d-%b-%Y %H:%M:%S",
            )
            for fmt in formats:
                try:
                    parsed = datetime.datetime.strptime(cleaned, fmt)
                    break
                except ValueError:
                    continue

        if parsed is None:
            raise ValueError("Could not parse HYCOM time origin %r" % text)
        if parsed.tzinfo is not None:
            parsed = parsed.astimezone(datetime.timezone.utc).replace(tzinfo=None)
        return parsed

    @classmethod
    def _decode_time(
        cls, variable: Any, values: np.ndarray
    ) -> Tuple[CurrentTimeList, datetime.datetime, float]:
        """Decode CF time units and return times, origin, seconds-per-unit."""

        attrs = cls._attributes(variable)
        units = str(attrs.get("units", "")).strip()
        match = re.match(
            r"^\s*(seconds?|minutes?|hours?|days?)\s+since\s+(.+?)\s*$",
            units,
            flags=re.IGNORECASE,
        )
        if match:
            unit_name = match.group(1).lower()
            origin = cls._parse_datetime(match.group(2))
            seconds_per_unit = {
                "second": 1.0,
                "seconds": 1.0,
                "minute": 60.0,
                "minutes": 60.0,
                "hour": 3600.0,
                "hours": 3600.0,
                "day": 86400.0,
                "days": 86400.0,
            }[unit_name]
        else:
            # Compatibility with the retired NCEI datasets, which exposed a
            # non-CF ``time_origin`` attribute and values in seconds.
            time_origin = attrs.get("time_origin")
            if time_origin is None:
                raise ValueError(
                    "HYCOM time variable has neither CF 'units since ...' "
                    "nor a 'time_origin' attribute"
                )
            origin = cls._parse_datetime(str(time_origin))
            seconds_per_unit = 1.0

        flat_values = np.asarray(values, dtype=float).reshape(-1)
        times = cls.time_list_class(
            [
                origin + datetime.timedelta(seconds=float(value) * seconds_per_unit)
                for value in flat_values
            ]
        )
        return times, origin, seconds_per_unit

    @staticmethod
    def _normalise_longitudes(longitudes: np.ndarray) -> np.ndarray:
        # Convert 0..360 to -180..180 without changing the source-grid order.
        result = ((np.asarray(longitudes, dtype=float) + 180.0) % 360.0) - 180.0
        # Avoid representing a source 180-degree point as -180 due only to
        # modulo arithmetic; either representation is equivalent spatially.
        return result

    @staticmethod
    def _nearest_index(values: Sequence[float], target: float) -> int:
        array = np.asarray(values, dtype=float)
        if not np.isfinite(array).any():
            raise ValueError("Cannot find a coordinate in an all-NaN axis")
        return int(np.nanargmin(np.abs(array - float(target))))

    @classmethod
    def _source_metadata(
        cls,
        dataset: Any,
        source_url: str,
        unit_variable: str,
        force: bool = False,
    ) -> Dict[str, Any]:
        if force:
            HYCOMDownloader._metadata_cache.pop(source_url, None)
        if source_url in HYCOMDownloader._metadata_cache:
            return HYCOMDownloader._metadata_cache[source_url]

        time_name, raw_time = cls._read_axis(dataset, "time")
        lat_name, raw_lats = cls._read_axis(dataset, "lat")
        lon_name, raw_lons = cls._read_axis(dataset, "lon")
        depth_name, raw_depths = cls._read_axis(dataset, "depth")

        time, ref_date, seconds_per_time_unit = cls._decode_time(
            dataset[time_name], raw_time
        )
        variable_name = cls._find_variable(dataset, (unit_variable,))
        units = str(cls._attributes(dataset[variable_name]).get("units", "m/s"))

        metadata = {
            "time_name": time_name,
            "lat_name": lat_name,
            "lon_name": lon_name,
            "depth_name": depth_name,
            "time": time,
            "raw_time": raw_time,
            "ref_date": ref_date,
            "seconds_per_time_unit": seconds_per_time_unit,
            "lats": np.asarray(raw_lats, dtype=float),
            "lons": cls._normalise_longitudes(raw_lons),
            # HyperNav historically represents ocean depth as negative.
            "depths": -np.asarray(raw_depths, dtype=float),
            "units": units,
        }
        HYCOMDownloader._metadata_cache[source_url] = metadata
        return metadata

    @classmethod
    def get_dimensions(
        cls,
        urlon: float,
        lllon: float,
        urlat: float,
        lllat: float,
        max_depth: float,
        dataset: Any,
    ):
        """Return regional axes and source indices using current CF names."""

        source_url = cls._resolve_source(cls.ID)
        metadata = cls._source_metadata(dataset, source_url, unit_variable="water_u")

        lower_lon_idx = cls._nearest_index(metadata["lons"], lllon)
        higher_lon_idx = cls._nearest_index(metadata["lons"], urlon)
        lower_lat_idx = cls._nearest_index(metadata["lats"], lllat)
        higher_lat_idx = cls._nearest_index(metadata["lats"], urlat)
        depth_idx = cls._nearest_index(metadata["depths"], max_depth)

        if lower_lon_idx > higher_lon_idx:
            raise ValueError(
                "This HYCOM regional reader does not support a bounding box "
                "that crosses the dateline"
            )
        if lower_lat_idx > higher_lat_idx:
            lower_lat_idx, higher_lat_idx = higher_lat_idx, lower_lat_idx

        # Slices use exclusive upper bounds.  Include the grid point nearest
        # the requested upper/right/deep boundaries.
        higher_lon_idx += 1
        higher_lat_idx += 1
        depth_stop = depth_idx + 1

        time = cls.time_list_class(list(metadata["time"]))
        lats = LatList(metadata["lats"][lower_lat_idx:higher_lat_idx].tolist())
        lons = LonList(metadata["lons"][lower_lon_idx:higher_lon_idx].tolist())
        depth = DepthList(metadata["depths"][:depth_stop].tolist())

        return (
            time,
            lats,
            lons,
            depth,
            lower_lon_idx,
            higher_lon_idx,
            lower_lat_idx,
            higher_lat_idx,
            metadata["units"],
            metadata["ref_date"],
        )

    @classmethod
    def initialize(cls, force: bool = False):
        """Connect this regional class to the current global HYCOM source."""

        if cls.__dict__.get("_initialized", False) and not force:
            return cls

        dataset = cls.get_dataset(cls.ID, force=force)
        dimensions = cls.get_dimensions(
            cls.urlon,
            cls.lllon,
            cls.urlat,
            cls.lllat,
            cls.max_depth,
            dataset,
        )
        (
            cls.dataset_time,
            cls.lats,
            cls.lons,
            cls.depths,
            cls.lllon_idx,
            cls.urlon_idx,
            cls.lllat_idx,
            cls.urlat_idx,
            cls.units,
            cls.ref_date,
        ) = dimensions
        cls.dataset = dataset
        cls._initialized = True
        return cls

    @classmethod
    def refresh(cls):
        """Drop velocity/tracer caches and reconnect to the rolling FMRC."""

        for source_url in {cls.velocity_url, cls.tracer_url}:
            HYCOMDownloader._dataset_cache.pop(source_url, None)
            HYCOMDownloader._metadata_cache.pop(source_url, None)
        cls._initialized = False
        return cls.initialize(force=True)

    @classmethod
    def source_time_range(
        cls, refresh_source: bool = False
    ) -> Tuple[datetime.datetime, datetime.datetime]:
        if refresh_source:
            cls.refresh()
        else:
            cls.initialize()
        return min(cls.dataset_time), max(cls.dataset_time)

    @classmethod
    def _read_grid(
        cls,
        dataset: Any,
        variable_name: str,
        index: Tuple[slice, slice, slice, slice],
    ) -> np.ndarray:
        result = dataset[variable_name][index]
        array = cls._unwrap_data(result, variable_name)
        return np.asanyarray(array)

    @staticmethod
    def _clean_numeric_array(array: np.ndarray) -> np.ndarray:
        data = np.ma.asarray(array, dtype=float)
        data = np.ma.masked_where(np.abs(data) > 1.0e10, data)
        return data.filled(np.nan)

    @classmethod
    def _chunk_slices(cls, chunk_size: int = 100):
        """Return ``(chunk_number, start, stop)`` slices for the source time axis."""

        chunk_size = int(chunk_size)
        if chunk_size < 1:
            raise ValueError("chunk_size must be a positive integer")
        n_times = len(cls.dataset_time)
        if n_times == 0:
            return []
        starts = list(range(0, n_times, chunk_size))
        return [
            (chunk_number, start, min(start + chunk_size, n_times))
            for chunk_number, start in enumerate(starts)
        ]

    @classmethod
    def _time_seconds(cls, start: int, stop: int) -> np.ndarray:
        return np.asarray(
            [
                (value - cls.ref_date).total_seconds()
                for value in cls.dataset_time[start:stop]
            ],
            dtype=float,
        )

    @staticmethod
    def _atomic_pickle(filename: str, payload: Mapping[str, Any]):
        """Write a pickle atomically so an interrupted transfer is never accepted."""

        folder = os.path.dirname(filename)
        if folder and not os.path.exists(folder):
            os.makedirs(folder)
        temporary = "%s.tmp-%d" % (filename, os.getpid())
        try:
            with open(temporary, "wb") as handle:
                pickle.dump(payload, handle, protocol=pickle.HIGHEST_PROTOCOL)
            os.replace(temporary, filename)
        finally:
            if os.path.exists(temporary):
                os.remove(temporary)

    @staticmethod
    def _utc_now_iso() -> str:
        return (
            datetime.datetime.now(datetime.timezone.utc)
            .isoformat(timespec="seconds")
            .replace("+00:00", "Z")
        )

    @staticmethod
    def _manifest_time(value: Optional[datetime.datetime]) -> Optional[str]:
        if value is None:
            return None
        if value.tzinfo is None:
            value = value.replace(tzinfo=datetime.timezone.utc)
        else:
            value = value.astimezone(datetime.timezone.utc)
        return value.isoformat(timespec="seconds").replace("+00:00", "Z")

    @classmethod
    def manifest_filename(cls) -> str:
        """Return the JSON manifest stored beside the numeric chunk files."""

        folder = os.path.dirname(cls.make_k_filename(0))
        if folder and not os.path.exists(folder):
            os.makedirs(folder)
        return os.path.join(folder, "manifest.json")

    @staticmethod
    def _atomic_json(filename: str, payload: Mapping[str, Any]):
        folder = os.path.dirname(filename)
        if folder and not os.path.exists(folder):
            os.makedirs(folder)
        temporary = "%s.tmp-%d" % (filename, os.getpid())
        try:
            with open(temporary, "w", encoding="utf-8") as handle:
                json.dump(payload, handle, indent=2, sort_keys=True, allow_nan=False)
                handle.write("\n")
            os.replace(temporary, filename)
        finally:
            if os.path.exists(temporary):
                os.remove(temporary)

    @classmethod
    def _new_download_manifest(cls) -> Dict[str, Any]:
        now = cls._utc_now_iso()
        return {
            "schema_version": 1,
            "dataset_description": cls.dataset_description,
            "model_name": cls.model_name,
            "region": cls.location,
            "created_at": now,
            "updated_at": now,
            "chunks": {},
        }

    @classmethod
    def _load_download_manifest(cls) -> Dict[str, Any]:
        """Read the manifest, recovering cleanly from an invalid JSON file."""

        filename = cls.manifest_filename()
        if not os.path.isfile(filename):
            return cls._new_download_manifest()
        try:
            with open(filename, "r", encoding="utf-8") as handle:
                manifest = json.load(handle)
            if not isinstance(manifest, dict):
                raise ValueError("manifest root is not a JSON object")
            if not isinstance(manifest.get("chunks", {}), dict):
                raise ValueError("manifest 'chunks' member is not a JSON object")
            manifest.setdefault("chunks", {})
            manifest.setdefault("created_at", cls._utc_now_iso())
            return manifest
        except (OSError, ValueError, TypeError, json.JSONDecodeError) as exc:
            timestamp = datetime.datetime.now(datetime.timezone.utc).strftime(
                "%Y%m%dT%H%M%SZ"
            )
            backup = "%s.corrupt-%s" % (filename, timestamp)
            try:
                os.replace(filename, backup)
                detail = " The invalid copy was moved to %s." % backup
            except OSError:
                detail = " The invalid copy could not be moved."
            warnings.warn(
                "Could not read the HYCOM download manifest: %s.%s "
                "A new manifest will be created." % (exc, detail),
                RuntimeWarning,
            )
            return cls._new_download_manifest()

    @classmethod
    def _source_fingerprint(cls) -> str:
        """Identify the exact rolling time window represented by this run."""

        source_url = cls._resolve_source(cls.ID)
        if not cls.dataset_time:
            return "%s|0" % source_url
        return "%s|%d|%s|%s" % (
            source_url,
            len(cls.dataset_time),
            cls._manifest_time(min(cls.dataset_time)),
            cls._manifest_time(max(cls.dataset_time)),
        )

    @classmethod
    def _expected_chunk_entry(
        cls,
        k: int,
        start: int,
        stop: int,
        source_fingerprint: Optional[str] = None,
    ) -> Dict[str, Any]:
        filename = cls.make_k_filename(k)
        times = list(cls.dataset_time[start:stop])
        return {
            "chunk": int(k),
            "start_index": int(start),
            "stop_index": int(stop),
            "record_count": int(stop - start),
            "start_time": cls._manifest_time(times[0]) if times else None,
            "end_time": cls._manifest_time(times[-1]) if times else None,
            "expected_shape": [
                int(stop - start),
                int(len(cls.depths)),
                int(len(cls.lats)),
                int(len(cls.lons)),
            ],
            "file": os.path.basename(filename),
            # The loader always reconstructs the path with make_k_filename();
            # this value is descriptive and never trusted as the source of truth.
            "path": os.path.abspath(filename),
            "source_fingerprint": source_fingerprint or cls._source_fingerprint(),
        }

    @staticmethod
    def _manifest_entry_matches_expected(
        entry: Mapping[str, Any], expected: Mapping[str, Any]
    ) -> bool:
        fields = (
            "chunk",
            "start_index",
            "stop_index",
            "record_count",
            "start_time",
            "end_time",
            "expected_shape",
            "file",
            "source_fingerprint",
        )
        return all(entry.get(field) == expected.get(field) for field in fields)

    @classmethod
    def _manifest_status_summary(cls, manifest: Mapping[str, Any]) -> Dict[str, Any]:
        chunks = manifest.get("chunks", {})
        expected_keys = manifest.get("expected_chunks", list(chunks))
        counts: Dict[str, int] = {}
        for key in expected_keys:
            status = str(chunks.get(str(key), {}).get("status", "pending"))
            counts[status] = counts.get(status, 0) + 1
        total = len(expected_keys)
        complete = counts.get("complete", 0)
        return {
            "expected": int(total),
            "complete": int(complete),
            "pending": int(counts.get("pending", 0)),
            "stale": int(counts.get("stale", 0)),
            "missing": int(counts.get("missing", 0)),
            "invalid": int(counts.get("invalid", 0)),
            "downloading": int(counts.get("downloading", 0)),
            "retrying": int(counts.get("retrying", 0)),
            "error": int(counts.get("error", 0)),
            "all_complete": bool(complete == total),
        }

    @classmethod
    def _write_download_manifest(cls, manifest: Dict[str, Any]):
        manifest["schema_version"] = 1
        manifest["dataset_description"] = cls.dataset_description
        manifest["model_name"] = cls.model_name
        manifest["region"] = cls.location
        manifest["updated_at"] = cls._utc_now_iso()
        manifest["summary"] = cls._manifest_status_summary(manifest)
        cls._atomic_json(cls.manifest_filename(), manifest)

    @classmethod
    def _prepare_download_manifest(
        cls,
        chunks: Sequence[Tuple[int, int, int]],
        chunk_size: int,
        prune_unexpected: bool,
    ) -> Dict[str, Any]:
        """Synchronize manifest metadata with the currently opened HYCOM source."""

        manifest = cls._load_download_manifest()
        source_fingerprint = cls._source_fingerprint()
        manifest.update(
            {
                "source_url": cls._resolve_source(cls.ID),
                "source_fingerprint": source_fingerprint,
                "source_time_start": (
                    cls._manifest_time(min(cls.dataset_time))
                    if cls.dataset_time
                    else None
                ),
                "source_time_end": (
                    cls._manifest_time(max(cls.dataset_time))
                    if cls.dataset_time
                    else None
                ),
                "source_record_count": int(len(cls.dataset_time)),
                "chunk_size": int(chunk_size),
            }
        )

        entries = manifest.setdefault("chunks", {})
        expected_keys = []
        for k, start, stop in chunks:
            key = str(k)
            expected_keys.append(key)
            expected = cls._expected_chunk_entry(
                k, start, stop, source_fingerprint=source_fingerprint
            )
            previous = entries.get(key)
            if isinstance(previous, Mapping) and cls._manifest_entry_matches_expected(
                previous, expected
            ):
                merged = dict(previous)
                merged.update(expected)
            else:
                merged = dict(expected)
                if isinstance(previous, Mapping):
                    merged["previous_status"] = previous.get("status")
                    merged["status"] = "stale"
                    merged["status_reason"] = (
                        "The source window or expected chunk metadata changed."
                    )
                else:
                    merged["status"] = "pending"
            entries[key] = merged

        if prune_unexpected:
            for key in list(entries):
                if key not in expected_keys:
                    entries.pop(key, None)

        manifest["expected_chunks"] = expected_keys
        manifest["expected_chunk_count"] = len(expected_keys)
        cls._write_download_manifest(manifest)
        return manifest

    @classmethod
    def _validate_chunk_file(
        cls,
        filename: str,
        start: int,
        stop: int,
        expected_shape: Sequence[int],
    ) -> Tuple[bool, str, Dict[str, Any]]:
        """Verify that a local chunk can be read and matches its expected data."""

        if not os.path.isfile(filename):
            return False, "file is missing", {}
        try:
            file_size = int(os.path.getsize(filename))
        except OSError as exc:
            return False, "could not stat file: %s" % exc, {}
        if file_size <= 0:
            return False, "file is empty", {"bytes": file_size}

        try:
            with open(filename, "rb") as handle:
                payload = pickle.load(handle)
        except Exception as exc:
            return (
                False,
                "pickle could not be read: %s: %s" % (type(exc).__name__, exc),
                {"bytes": file_size},
            )
        if not isinstance(payload, Mapping):
            return False, "pickle payload is not a mapping", {"bytes": file_size}

        missing_keys = [key for key in ("u", "v", "time") if key not in payload]
        if missing_keys:
            return (
                False,
                "pickle is missing keys: %s" % ", ".join(missing_keys),
                {"bytes": file_size},
            )

        expected_shape_tuple = tuple(int(value) for value in expected_shape)
        try:
            u = np.asanyarray(payload["u"])
            v = np.asanyarray(payload["v"])
            time_values = np.asarray(payload["time"], dtype=float).reshape(-1)
        except Exception as exc:
            return (
                False,
                "payload arrays could not be decoded: %s: %s"
                % (type(exc).__name__, exc),
                {"bytes": file_size},
            )

        if u.shape != expected_shape_tuple or v.shape != expected_shape_tuple:
            return (
                False,
                "array shapes are u=%r, v=%r; expected=%r"
                % (u.shape, v.shape, expected_shape_tuple),
                {"bytes": file_size},
            )
        if u.shape != v.shape:
            return False, "u and v shapes differ", {"bytes": file_size}

        expected_time = cls._time_seconds(start, stop)
        if time_values.shape != expected_time.shape:
            return (
                False,
                "time shape is %r; expected %r"
                % (time_values.shape, expected_time.shape),
                {"bytes": file_size},
            )
        if not np.all(np.isfinite(time_values)):
            return False, "time coordinate contains non-finite values", {
                "bytes": file_size
            }
        if not np.allclose(time_values, expected_time, rtol=0.0, atol=1.0e-6):
            return (
                False,
                "time coordinate does not match the current source window",
                {"bytes": file_size},
            )

        modified = datetime.datetime.fromtimestamp(
            os.path.getmtime(filename), tz=datetime.timezone.utc
        )
        return (
            True,
            "",
            {
                "bytes": file_size,
                "modified_at": cls._manifest_time(modified),
                "actual_shape": list(u.shape),
                "time_records": int(time_values.size),
            },
        )

    @classmethod
    def _download_chunks(
        cls,
        chunks: Iterable[Tuple[int, int, int]],
        overwrite: bool,
        max_attempts: int,
        manifest: Dict[str, Any],
        verify_existing: bool = True,
        operation: str = "download",
    ) -> Tuple[Dict[str, Any], Dict[str, Any]]:
        """Download selected chunks and persist state after every attempt.

        Existing files are skipped only after they have been validated.  Any
        manifest entry with an interrupted, stale, invalid, or error state is
        downloaded again even if a file with the expected name is present.
        """

        cls.initialize()
        max_attempts = int(max_attempts)
        if max_attempts < 1:
            raise ValueError("max_attempts must be at least one")

        retry_errors = (
            HTTPError,
            SocketTimeout,
            OSError,
            RuntimeError,
            ValueError,
            RequestException,
            DapError,
        )
        failed_states = {
            "error",
            "failed",
            "invalid",
            "missing",
            "retrying",
            "downloading",
            "stale",
        }
        chunks = list(chunks)
        run_started = cls._utc_now_iso()
        summary: Dict[str, Any] = {
            "operation": operation,
            "manifest": cls.manifest_filename(),
            "started_at": run_started,
            "downloaded": [],
            "redownloaded": [],
            "skipped": [],
            "invalid": {},
            "failed": {},
        }

        manifest["last_run"] = {
            "operation": operation,
            "started_at": run_started,
            "requested_chunks": [int(k) for k, _, _ in chunks],
            "overwrite": bool(overwrite),
            "verify_existing": bool(verify_existing),
            "max_attempts": int(max_attempts),
            "status": "running",
        }
        cls._write_download_manifest(manifest)

        for k, start, stop in chunks:
            key = str(k)
            filename = cls.make_k_filename(k)
            expected = cls._expected_chunk_entry(
                k,
                start,
                stop,
                source_fingerprint=str(manifest.get("source_fingerprint", "")),
            )
            previous = manifest.setdefault("chunks", {}).get(key, {})
            if not isinstance(previous, Mapping):
                previous = {}
            entry = dict(previous)
            entry.update(expected)
            previous_status = str(entry.get("status", "pending"))
            file_existed_before = os.path.isfile(filename)

            # A failed/interrupted/stale entry is deliberately retried.  Other
            # files may be skipped, but only if they validate against this exact
            # source time window and expected array shape.
            force_redownload = bool(overwrite or previous_status in failed_states)
            if file_existed_before and not force_redownload:
                if verify_existing:
                    valid, reason, details = cls._validate_chunk_file(
                        filename,
                        start,
                        stop,
                        expected["expected_shape"],
                    )
                else:
                    valid = previous_status == "complete"
                    reason = "manifest does not mark the file complete"
                    details = {
                        "bytes": int(os.path.getsize(filename))
                        if os.path.isfile(filename)
                        else 0
                    }

                if valid:
                    entry.update(details)
                    entry.update(
                        {
                            "status": "complete",
                            "last_action": "skipped_existing_valid_file",
                            "validated_at": cls._utc_now_iso(),
                            "last_error": None,
                            "last_error_type": None,
                        }
                    )
                    manifest["chunks"][key] = entry
                    summary["skipped"].append(int(k))
                    print("[skip] HYCOM chunk %d is already complete: %s" % (k, filename))
                    cls._write_download_manifest(manifest)
                    continue

                entry.update(
                    {
                        "status": "invalid",
                        "status_reason": reason,
                        "invalidated_at": cls._utc_now_iso(),
                    }
                )
                entry.update(details)
                manifest["chunks"][key] = entry
                summary["invalid"][key] = reason
                print("[repair] HYCOM chunk %d is invalid (%s); downloading again" % (k, reason))
                cls._write_download_manifest(manifest)
                force_redownload = True
            elif not file_existed_before and previous_status == "complete":
                entry.update(
                    {
                        "status": "missing",
                        "status_reason": "manifest marked complete but file is missing",
                        "invalidated_at": cls._utc_now_iso(),
                    }
                )
                manifest["chunks"][key] = entry
                summary["invalid"][key] = entry["status_reason"]
                print("[repair] HYCOM chunk %d is missing; downloading again" % k)
                cls._write_download_manifest(manifest)
                force_redownload = True

            success = False
            last_exception: Optional[Exception] = None
            for attempt in range(1, max_attempts + 1):
                attempt_started = cls._utc_now_iso()
                entry.update(
                    {
                        "status": "downloading",
                        "last_action": "download",
                        "attempt_in_run": int(attempt),
                        "total_attempts": int(entry.get("total_attempts", 0)) + 1,
                        "last_attempt_started_at": attempt_started,
                    }
                )
                manifest["chunks"][key] = entry
                cls._write_download_manifest(manifest)

                try:
                    index = (
                        slice(start, stop),
                        slice(0, len(cls.depths)),
                        slice(cls.lllat_idx, cls.urlat_idx),
                        slice(cls.lllon_idx, cls.urlon_idx),
                    )
                    u = cls._clean_numeric_array(
                        cls._read_grid(cls.dataset, "water_u", index)
                    )
                    v = cls._clean_numeric_array(
                        cls._read_grid(cls.dataset, "water_v", index)
                    )
                    time_values = cls._time_seconds(start, stop)
                    expected_shape = tuple(expected["expected_shape"])
                    if u.shape != expected_shape or v.shape != expected_shape:
                        raise ValueError(
                            "Unexpected HYCOM subset shape: u=%r, v=%r, expected=%r"
                            % (u.shape, v.shape, expected_shape)
                        )

                    cls._atomic_pickle(
                        filename,
                        {"u": u, "v": v, "time": time_values},
                    )
                    valid, reason, details = cls._validate_chunk_file(
                        filename,
                        start,
                        stop,
                        expected["expected_shape"],
                    )
                    if not valid:
                        raise ValueError(
                            "newly written chunk failed validation: %s" % reason
                        )

                    entry.update(details)
                    entry.update(
                        {
                            "status": "complete",
                            "downloaded_at": cls._utc_now_iso(),
                            "validated_at": cls._utc_now_iso(),
                            "attempts_last_run": int(attempt),
                            "last_error": None,
                            "last_error_type": None,
                            "status_reason": None,
                        }
                    )
                    manifest["chunks"][key] = entry
                    cls._write_download_manifest(manifest)
                    summary["downloaded"].append(int(k))
                    if file_existed_before or force_redownload:
                        summary["redownloaded"].append(int(k))
                    print("[ok] HYCOM chunk %d downloaded and validated: %s" % (k, filename))
                    success = True
                    break
                except Exception as exc:  # persist every failure in the manifest
                    last_exception = exc
                    retryable = isinstance(exc, retry_errors)
                    will_retry = bool(retryable and attempt < max_attempts)
                    entry.update(
                        {
                            "status": "retrying" if will_retry else "error",
                            "last_error": str(exc),
                            "last_error_type": type(exc).__name__,
                            "last_error_at": cls._utc_now_iso(),
                            "attempts_last_run": int(attempt),
                        }
                    )
                    manifest["chunks"][key] = entry
                    cls._write_download_manifest(manifest)

                    if not will_retry:
                        break
                    warnings.warn(
                        "HYCOM chunk %d failed on attempt %d; reconnecting: %s"
                        % (k, attempt, exc),
                        RuntimeWarning,
                    )
                    try:
                        # Keep the original axis/time snapshot.  The manifest
                        # and payload-time validation prevent a silently shifted
                        # rolling FMRC window from being accepted as complete.
                        cls.dataset = cls.get_dataset(cls.ID, force=True)
                    except Exception as reconnect_exc:
                        warnings.warn(
                            "HYCOM reconnect after chunk %d failed: %s"
                            % (k, reconnect_exc),
                            RuntimeWarning,
                        )

            if not success:
                message = (
                    "%s: %s" % (type(last_exception).__name__, last_exception)
                    if last_exception is not None
                    else "unknown download error"
                )
                summary["failed"][key] = message
                print("[failed] HYCOM chunk %d: %s" % (k, message))

        completed_at = cls._utc_now_iso()
        manifest["last_run"].update(
            {
                "completed_at": completed_at,
                "status": "error" if summary["failed"] else "complete",
                "downloaded": summary["downloaded"],
                "redownloaded": summary["redownloaded"],
                "skipped": summary["skipped"],
                "failed": summary["failed"],
            }
        )
        cls._write_download_manifest(manifest)
        summary["completed_at"] = completed_at
        summary["manifest_summary"] = dict(manifest.get("summary", {}))
        return summary, manifest

    @classmethod
    def download_and_save(
        cls,
        overwrite: bool = False,
        chunk_size: int = 100,
        max_attempts: int = 3,
        refresh_source: bool = True,
        verify_existing: bool = True,
        raise_on_error: bool = True,
    ) -> Dict[str, Any]:
        """Download every current HYCOM chunk using a resumable manifest.

        The default behavior is safe to rerun:

        * valid files already on disk are skipped;
        * missing, stale, interrupted, corrupt, or previously failed chunks are
          downloaded again;
        * each attempt and final validation result is written atomically to
          ``manifest.json`` beside the chunk files.

        ``overwrite=True`` forces every current chunk to be downloaded again.
        The method attempts all chunks before raising a summarized error, so a
        failure in one chunk does not prevent later chunks from completing.
        """

        if refresh_source:
            cls.refresh()
        else:
            cls.initialize()
        chunks = cls._chunk_slices(chunk_size=chunk_size)
        manifest = cls._prepare_download_manifest(
            chunks,
            chunk_size=int(chunk_size),
            prune_unexpected=True,
        )
        summary, manifest = cls._download_chunks(
            chunks,
            overwrite=bool(overwrite),
            max_attempts=max_attempts,
            manifest=manifest,
            verify_existing=verify_existing,
            operation="download_all",
        )

        # Remove numeric chunks left by a previously longer rolling window.
        folder = os.path.dirname(cls.make_k_filename(0))
        valid = {str(k) for k, _, _ in chunks}
        removed = []
        if os.path.isdir(folder):
            for basename in os.listdir(folder):
                if basename.isdigit() and basename not in valid:
                    path = os.path.join(folder, basename)
                    os.remove(path)
                    removed.append(path)
        summary["removed_stale_files"] = removed
        manifest["last_run"]["removed_stale_files"] = removed
        cls._write_download_manifest(manifest)
        summary["manifest_summary"] = dict(manifest.get("summary", {}))

        if summary["failed"] and raise_on_error:
            raise RuntimeError(
                "%d HYCOM chunk(s) failed after all retries. See %s for details: %s"
                % (
                    len(summary["failed"]),
                    cls.manifest_filename(),
                    ", ".join(sorted(summary["failed"])),
                )
            )
        return summary

    @classmethod
    def download_recent(
        cls,
        recent_chunks: int = 2,
        chunk_size: int = 100,
        max_attempts: int = 3,
        refresh_source: bool = True,
        overwrite: bool = True,
        verify_existing: bool = True,
        raise_on_error: bool = True,
    ) -> Dict[str, Any]:
        """Download or refresh the newest chunks while updating the manifest."""

        if refresh_source:
            cls.refresh()
        else:
            cls.initialize()
        all_chunks = cls._chunk_slices(chunk_size=chunk_size)
        recent_chunks = max(1, int(recent_chunks))
        selected = all_chunks[max(0, len(all_chunks) - recent_chunks) :]
        manifest = cls._prepare_download_manifest(
            all_chunks,
            chunk_size=int(chunk_size),
            prune_unexpected=False,
        )
        summary, _ = cls._download_chunks(
            selected,
            overwrite=bool(overwrite),
            max_attempts=max_attempts,
            manifest=manifest,
            verify_existing=verify_existing,
            operation="download_recent",
        )
        if summary["failed"] and raise_on_error:
            raise RuntimeError(
                "%d recent HYCOM chunk(s) failed after all retries. See %s: %s"
                % (
                    len(summary["failed"]),
                    cls.manifest_filename(),
                    ", ".join(sorted(summary["failed"])),
                )
            )
        return summary

    @classmethod
    def audit_downloads(
        cls,
        chunk_size: int = 100,
        refresh_source: bool = False,
    ) -> Dict[str, Any]:
        """Validate all expected files and update their manifest statuses.

        This method performs no network data downloads.  A subsequent
        ``download_and_save(overwrite=False)`` repairs every entry marked
        ``missing``, ``invalid``, ``stale``, or ``error``.
        """

        if refresh_source:
            cls.refresh()
        else:
            cls.initialize()
        chunks = cls._chunk_slices(chunk_size=chunk_size)
        manifest = cls._prepare_download_manifest(
            chunks,
            chunk_size=int(chunk_size),
            prune_unexpected=True,
        )
        result: Dict[str, Any] = {
            "manifest": cls.manifest_filename(),
            "complete": [],
            "missing": [],
            "invalid": {},
        }
        for k, start, stop in chunks:
            key = str(k)
            filename = cls.make_k_filename(k)
            expected = cls._expected_chunk_entry(
                k,
                start,
                stop,
                source_fingerprint=str(manifest.get("source_fingerprint", "")),
            )
            entry = dict(manifest["chunks"].get(key, {}))
            entry.update(expected)
            valid, reason, details = cls._validate_chunk_file(
                filename,
                start,
                stop,
                expected["expected_shape"],
            )
            if valid:
                entry.update(details)
                entry.update(
                    {
                        "status": "complete",
                        "last_action": "audit",
                        "validated_at": cls._utc_now_iso(),
                        "last_error": None,
                        "last_error_type": None,
                    }
                )
                result["complete"].append(int(k))
            else:
                status = "missing" if reason == "file is missing" else "invalid"
                entry.update(
                    {
                        "status": status,
                        "status_reason": reason,
                        "invalidated_at": cls._utc_now_iso(),
                        "last_action": "audit",
                    }
                )
                if status == "missing":
                    result["missing"].append(int(k))
                else:
                    result["invalid"][key] = reason
            manifest["chunks"][key] = entry

        manifest["last_audit"] = {
            "completed_at": cls._utc_now_iso(),
            "complete": result["complete"],
            "missing": result["missing"],
            "invalid": result["invalid"],
        }
        cls._write_download_manifest(manifest)
        result["manifest_summary"] = dict(manifest.get("summary", {}))
        return result

    @staticmethod
    def _normalise_datetime(value: datetime.datetime) -> datetime.datetime:
        if value.tzinfo is None:
            return value
        return value.astimezone(datetime.timezone.utc).replace(tzinfo=None)

    @staticmethod
    def _normalise_spatial_stride(spatial_stride: Any) -> Tuple[int, int]:
        """Return ``(latitude_stride, longitude_stride)`` for a subset."""

        if isinstance(spatial_stride, (int, np.integer)):
            lat_stride = lon_stride = int(spatial_stride)
        else:
            try:
                lat_stride, lon_stride = spatial_stride
            except (TypeError, ValueError) as exc:
                raise ValueError(
                    "spatial_stride must be an integer or a two-item "
                    "(latitude, longitude) sequence"
                ) from exc
            lat_stride = int(lat_stride)
            lon_stride = int(lon_stride)
        if lat_stride < 1 or lon_stride < 1:
            raise ValueError("spatial strides must be positive integers")
        return lat_stride, lon_stride

    @classmethod
    def _coordinate_window(
        cls,
        values: Sequence[float],
        lower: float,
        upper: float,
        axis_name: str,
    ) -> Tuple[int, int]:
        """Return an inclusive-nearest coordinate window as slice bounds."""

        lower = float(lower)
        upper = float(upper)
        if lower > upper:
            if axis_name == "longitude":
                raise ValueError(
                    "Longitude bounds must increase west-to-east; bounding "
                    "boxes that cross the dateline are not supported"
                )
            raise ValueError("%s lower bound exceeds upper bound" % axis_name)

        start = cls._nearest_index(values, lower)
        stop = cls._nearest_index(values, upper)
        if start > stop:
            # This can happen only across the discontinuity introduced when
            # HYCOM's 0..360 longitude axis is represented as -180..180.
            if axis_name == "longitude":
                raise ValueError(
                    "The requested longitude window crosses HYCOM's dateline "
                    "discontinuity and cannot be represented by one slice"
                )
            start, stop = stop, start
        return start, stop + 1

    @staticmethod
    def _reshape_grid_subset(
        array: np.ndarray,
        expected_shape: Tuple[int, int, int, int],
        variable_name: str,
    ) -> np.ndarray:
        """Restore singleton dimensions that an OPeNDAP client may squeeze."""

        array = np.asarray(array, dtype=float)
        expected_size = int(np.prod(expected_shape))
        if array.size != expected_size:
            raise ValueError(
                "Unexpected %s subset shape %r; expected %r"
                % (variable_name, array.shape, expected_shape)
            )
        return array.reshape(expected_shape)

    @staticmethod
    def _teos10_sound_speed(
        salinity: np.ndarray,
        temperature: np.ndarray,
        depths: Sequence[float],
        latitudes: Sequence[float],
        longitudes: Sequence[float],
        temperature_kind: str = "in_situ",
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """Convert HYCOM tracers to TEOS-10 sound speed.

        HYCOM's published ``water_temp`` variable has the CF standard name
        ``sea_water_temperature``, so ``in_situ`` is the default.  Set
        ``temperature_kind='potential'`` only when using a custom HYCOM source
        whose temperature is potential temperature referenced to 0 dbar, or
        ``temperature_kind='conservative'`` when it already contains TEOS-10
        Conservative Temperature.
        """

        practical_salinity = np.asarray(salinity, dtype=float)
        water_temperature = np.asarray(temperature, dtype=float)
        if practical_salinity.shape != water_temperature.shape:
            raise ValueError("salinity and temperature arrays must have one shape")
        if practical_salinity.ndim != 3:
            raise ValueError(
                "salinity and temperature must have (depth, latitude, longitude) axes"
            )

        target_shape = practical_salinity.shape
        z = np.broadcast_to(
            np.asarray(depths, dtype=float)[:, np.newaxis, np.newaxis],
            target_shape,
        )
        lat = np.broadcast_to(
            np.asarray(latitudes, dtype=float)[np.newaxis, :, np.newaxis],
            target_shape,
        )
        lon = np.broadcast_to(
            np.asarray(longitudes, dtype=float)[np.newaxis, np.newaxis, :],
            target_shape,
        )

        import gsw

        pressure = gsw.p_from_z(z, lat)
        absolute_salinity = gsw.SA_from_SP(practical_salinity, pressure, lon, lat)

        kind = str(temperature_kind).strip().lower().replace("-", "_")
        if kind in {"in_situ", "insitu", "t"}:
            conservative_temperature = gsw.CT_from_t(
                absolute_salinity, water_temperature, pressure
            )
        elif kind in {"potential", "potential_temperature", "pt"}:
            conservative_temperature = gsw.CT_from_pt(
                absolute_salinity, water_temperature
            )
        elif kind in {"conservative", "conservative_temperature", "ct"}:
            conservative_temperature = water_temperature
        else:
            raise ValueError(
                "temperature_kind must be 'in_situ', 'potential', or 'conservative'"
            )

        sound_speed = gsw.sound_speed(
            absolute_salinity, conservative_temperature, pressure
        )
        invalid = ~np.isfinite(practical_salinity) | ~np.isfinite(water_temperature)
        sound_speed = np.array(sound_speed, dtype=float, copy=True)
        absolute_salinity = np.array(absolute_salinity, dtype=float, copy=True)
        conservative_temperature = np.array(
            conservative_temperature, dtype=float, copy=True
        )
        pressure = np.array(pressure, dtype=float, copy=True)
        sound_speed[invalid] = np.nan
        absolute_salinity[invalid] = np.nan
        conservative_temperature[invalid] = np.nan
        pressure[invalid] = np.nan
        return (
            sound_speed,
            absolute_salinity,
            conservative_temperature,
            pressure,
        )

    @classmethod
    def get_sound_speed_profiles(
        cls,
        when: Optional[datetime.datetime] = None,
        lllon: Optional[float] = None,
        urlon: Optional[float] = None,
        lllat: Optional[float] = None,
        urlat: Optional[float] = None,
        max_depth: Optional[float] = None,
        spatial_stride: Any = 1,
        temperature_kind: str = "in_situ",
        drop_empty: bool = True,
        refresh_source: bool = False,
    ) -> Dict[str, Any]:
        """Return every vertical sound-speed profile in a regional box.

        Parameters
        ----------
        when
            Model time to use.  The nearest three-hourly field is selected;
            ``None`` selects the field nearest the current UTC time, clipped
            to the rolling tracer feed's available interval.
        lllon, urlon, lllat, urlat
            West, east, south, and north bounds in conventional -180..180
            longitude coordinates.  Omitted bounds use the regional class's
            configured box.
        max_depth
            Deepest profile level.  Positive or negative metres are accepted;
            the returned depth coordinate follows HyperNav's negative-down
            convention.  ``None`` uses ``cls.max_depth``.
        spatial_stride
            One positive integer for both horizontal axes, or
            ``(latitude_stride, longitude_stride)``.  A value of one retains
            all HYCOM grid points in the box.
        temperature_kind
            Interpretation of ``water_temp`` for TEOS-10 conversion.  The
            current ESPC-D-V02 CF field is in-situ temperature.
        refresh_source
            Reopen the rolling tracer endpoint and refresh its coordinate/time
            metadata before selecting the field.

        Returns
        -------
        dict
            ``sound_speed`` has shape ``(profile, depth)`` in m/s.  Parallel
            ``latitudes`` and ``longitudes`` arrays identify every profile.
        """

        tracer_dataset = cls.get_tracer_dataset(force=refresh_source)
        metadata = cls._source_metadata(
            tracer_dataset,
            cls.tracer_url,
            unit_variable="water_temp",
            force=refresh_source,
        )
        tracer_times = list(metadata["time"])
        if not tracer_times:
            raise ValueError("The HYCOM tracer source contains no time records")

        available_start, available_end = min(tracer_times), max(tracer_times)
        if when is None:
            current_utc = datetime.datetime.now(datetime.timezone.utc).replace(
                tzinfo=None
            )
            when = min(max(current_utc, available_start), available_end)
        else:
            when = cls._normalise_datetime(when)
        if when < available_start or when > available_end:
            raise ValueError(
                "Requested sound-speed time %s is outside the rolling "
                "ESPC-D-V02 tracer window %s to %s"
                % (when, available_start, available_end)
            )
        time_idx = min(
            range(len(tracer_times)),
            key=lambda index: abs(tracer_times[index] - when),
        )

        west = cls.lllon if lllon is None else float(lllon)
        east = cls.urlon if urlon is None else float(urlon)
        south = cls.lllat if lllat is None else float(lllat)
        north = cls.urlat if urlat is None else float(urlat)
        lon_start, lon_stop = cls._coordinate_window(
            metadata["lons"], west, east, "longitude"
        )
        lat_start, lat_stop = cls._coordinate_window(
            metadata["lats"], south, north, "latitude"
        )
        lat_stride, lon_stride = cls._normalise_spatial_stride(spatial_stride)

        target_depth = cls.max_depth if max_depth is None else -abs(float(max_depth))
        depth_stop = cls._nearest_index(metadata["depths"], target_depth) + 1
        selected_depths = np.asarray(metadata["depths"][:depth_stop], dtype=float)
        selected_lats = np.asarray(
            metadata["lats"][lat_start:lat_stop:lat_stride], dtype=float
        )
        selected_lons = np.asarray(
            metadata["lons"][lon_start:lon_stop:lon_stride], dtype=float
        )
        if selected_lats.size == 0 or selected_lons.size == 0:
            raise ValueError("The requested sound-speed region contains no grid cells")

        index = (
            slice(time_idx, time_idx + 1),
            slice(0, depth_stop),
            slice(lat_start, lat_stop, lat_stride),
            slice(lon_start, lon_stop, lon_stride),
        )
        expected_shape = (
            1,
            len(selected_depths),
            len(selected_lats),
            len(selected_lons),
        )
        salinity = cls._reshape_grid_subset(
            cls._clean_numeric_array(cls._read_grid(tracer_dataset, "salinity", index)),
            expected_shape,
            "salinity",
        )[0]
        temperature = cls._reshape_grid_subset(
            cls._clean_numeric_array(
                cls._read_grid(tracer_dataset, "water_temp", index)
            ),
            expected_shape,
            "water_temp",
        )[0]

        sound_speed, _, _, _ = cls._teos10_sound_speed(
            salinity,
            temperature,
            selected_depths,
            selected_lats,
            selected_lons,
            temperature_kind=temperature_kind,
        )

        lon_grid, lat_grid = np.meshgrid(selected_lons, selected_lats)
        profiles = sound_speed.transpose(1, 2, 0).reshape(-1, len(selected_depths))
        profile_lats = lat_grid.reshape(-1)
        profile_lons = lon_grid.reshape(-1)
        if drop_empty:
            valid = np.any(np.isfinite(profiles), axis=1)
            profiles = profiles[valid]
            profile_lats = profile_lats[valid]
            profile_lons = profile_lons[valid]
        if profiles.shape[0] == 0:
            raise ValueError(
                "No wet HYCOM sound-speed profiles were found in the requested box"
            )

        return {
            "sound_speed": profiles,
            "depths": selected_depths,
            "latitudes": profile_lats,
            "longitudes": profile_lons,
            "time": tracer_times[time_idx],
            "requested_time": when,
            "bounds": (west, east, south, north),
            "spatial_stride": (lat_stride, lon_stride),
            "temperature_kind": temperature_kind,
            "source_url": cls.tracer_url,
        }



def opendap_subset(
    session, ncss_url, *, variables, times, depths, grid, timeout, source_axes=None
):
    """Load only a contiguous, exact HYCOM velocity slab through OPeNDAP.

    Raw longitude conventions are matched geometrically without sorting the
    remote fields. The caller validates and normalizes the loaded result.
    The shared requests session uses bounded connection/read timeouts.
    """
    source_url = urlsplit(str(ncss_url))
    prefix = "/thredds/ncss/grid/"
    if (
        source_url.scheme != "https"
        or source_url.netloc != "ncss.hycom.org"
        or not source_url.path.startswith(prefix)
        or not source_url.path[len(prefix) :]
        or source_url.query
        or source_url.fragment
    ):
        raise ValueError("Unsupported HYCOM NCSS archive URL for OPeNDAP.")
    url = "dap2://tds.hycom.org/thredds/dodsC/" + source_url.path[len(prefix) :]
    if (
        isinstance(timeout, (bool, np.bool_))
        or not np.isfinite(timeout)
        or timeout <= 0
    ):
        raise ValueError("OPeNDAP timeout must be finite and positive.")
    variables = (variables,) if isinstance(variables, str) else tuple(variables)
    if (
        not variables
        or len(set(variables)) != len(variables)
        or any(not isinstance(name, str) or not name for name in variables)
    ):
        raise ValueError("OPeNDAP variables must be distinct nonempty names.")
    requested_times = pd.DatetimeIndex(pd.to_datetime(times, utc=True)).as_unit("ns")
    if requested_times.empty or requested_times.hasnans:
        raise ValueError("OPeNDAP needs nonmissing requested timestamps.")
    if set(grid) != {"lat", "lon"}:
        raise ValueError("OPeNDAP grid must supply exact lat and lon nodes.")
    wanted = {"time": requested_times.asi8, "depth": depths, **grid}
    with xr.open_dataset(
        url,
        engine="pydap",
        create_default_indexes=False,
        decode_times=False,
        drop_variables=["tau"],  # HYCOM's forecast-hour units are not CF dates.
        backend_kwargs={"session": session, "timeout": (20, float(timeout))},
    ) as source:
        dimensions = {"time", "depth", "lat", "lon"}
        if any(
            name not in source.data_vars or set(source[name].dims) != dimensions
            for name in variables
        ):
            raise ValueError("Requested HYCOM velocities must span time/depth/lat/lon.")
        indexers = {}
        for name, requested in wanted.items():
            if name not in source or source[name].dims != (name,):
                raise ValueError(
                    f"OPeNDAP {name} must be a one-dimensional source axis."
                )
            if source_axes is None:
                # Compatibility for direct callers without NCSS metadata.
                # Decode only this axis; the ordinary archive path below
                # already has verified native axes and never reads them whole.
                raw = np.asarray(xr.decode_cf(source[[name]].load())[name].values)
            else:
                if name not in source_axes:
                    raise ValueError(f"Source metadata lacks the OPeNDAP {name} axis.")
                raw = np.asarray(source_axes[name])
                if raw.ndim != 1 or len(raw) != source.sizes[name]:
                    raise ValueError(
                        f"OPeNDAP {name} shape differs from source metadata."
                    )
            if name == "time":
                if isinstance(source_axes, dict) and isinstance(
                    source_axes[name], pd.DatetimeIndex
                ):
                    raw = (
                        pd.to_datetime(source_axes[name], utc=True)
                        .tz_localize(None)
                        .to_numpy()
                    )
                if not np.issubdtype(raw.dtype, np.datetime64) or np.isnat(raw).any():
                    raise ValueError(
                        "OPeNDAP time must decode to nonmissing datetime64."
                    )
                values = pd.DatetimeIndex(raw).as_unit("ns").asi8
            else:
                if (
                    not np.issubdtype(raw.dtype, np.number)
                    or not np.isfinite(raw).all()
                ):
                    raise ValueError(
                        f"OPeNDAP {name} needs finite numeric source nodes."
                    )
                values = (raw + 180) % 360 - 180 if name == "lon" else raw
            requested = np.asarray(requested)
            if (
                requested.ndim != 1
                or not requested.size
                or not np.isfinite(requested).all()
                or np.any(np.diff(requested) <= 0)
            ):
                raise ValueError(
                    f"Requested OPeNDAP {name} nodes must increase uniquely."
                )
            order = np.argsort(values)
            ordered = values[order]
            if not len(ordered) or np.any(np.diff(ordered) <= 0):
                raise ValueError(
                    f"OPeNDAP {name} source nodes must be nonempty and unique."
                )
            insertion = np.searchsorted(ordered, requested)
            lower = np.clip(insertion - 1, 0, len(ordered) - 1)
            upper = np.clip(insertion, 0, len(ordered) - 1)
            closest = np.where(
                np.abs(ordered[lower] - requested)
                <= np.abs(ordered[upper] - requested),
                lower,
                upper,
            )
            matched = ordered[closest]
            exact = (
                matched == requested
                if name == "time"
                else np.isclose(matched, requested, rtol=0, atol=1e-6)
            )
            if not exact.all():
                raise ValueError(f"Requested OPeNDAP {name} nodes are unavailable.")
            indices = order[closest]
            stride = int(indices[1] - indices[0]) if len(indices) > 1 else 1
            if (
                stride <= 0
                or (len(indices) > 1 and not np.all(np.diff(indices) == stride))
                or (name != "time" and stride != 1)
            ):
                raise ValueError(
                    f"Non-contiguous OPeNDAP {name} selection is unsupported."
                )
            indexers[name] = slice(int(indices[0]), int(indices[-1]) + 1, stride)
        # Selecting variables and applying integer slabs before .load() makes
        # PyDAP send constrained regional DODS requests, never global fields.
        return xr.decode_cf(source[list(variables)].isel(indexers).load())
