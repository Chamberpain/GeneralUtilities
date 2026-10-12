"""Save source-native ROMS C-grid subsets without decoding or regridding fields.

Callers supply approved output and ComputePathHandler scratch directories. Imports
perform no I/O. Each archive has one immutable request identity; each finalized
NetCDF retains raw arrays/attributes and is checked before atomic publication.
"""
from __future__ import annotations

from contextlib import contextmanager
from dataclasses import dataclass
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import tempfile
import threading
import time
from urllib.parse import urljoin, urlsplit

import numpy as np
import pandas as pd
import xarray as xr

from .current_archive import NETCDF_LOCK

UCSC_URL = (
    "https://oceanmodeling.ucsc.edu/thredds/dodsC/"
    "ccsra_2016a_phys_agg_slevs/fmrc/"
    "CCSRA_2016a_Phys_ROMS_Sigma-level_Aggregation_best.ncd"
)
CORE_FIELDS = ("u", "v", "w", "zeta", "ubar", "vbar")
REQUIRED_FIELDS = (
    "u", "v", "lat_rho", "lon_rho", "h", "angle", "s_rho", "s_w",
    "Cs_r", "Cs_w", "hc",
)
GRID_FIELDS = (
    "s_rho", "s_w", "Cs_r", "Cs_w", "h", "hc", "angle", "pm", "pn",
    "Vtransform", "Vstretching", "theta_s", "theta_b", "Tcline", "spherical",
    "grid", *(f"{kind}_{grid}" for grid in ("rho", "u", "v", "psi")
              for kind in ("lat", "lon", "mask")),
)
_LOCKS: dict[str, threading.RLock] = {}
_LOCKS_GUARD = threading.Lock()


class ArchiveConflictError(ValueError):
    """An existing archive/file does not match the requested provenance."""


class MissingNativeTime(ValueError):
    """The requested time is absent; never substitute another source time."""


def _utc(value):
    result = pd.Timestamp(value)
    if pd.isna(result):
        raise ValueError("A finite UTC timestamp is required.")
    return result.tz_convert("UTC").tz_localize(None) if result.tz else result


def _json_value(value):
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, (Path, pd.Timestamp)):
        return str(value)
    raise TypeError(f"Cannot encode {type(value).__name__} as archive metadata.")


def _json_bytes(value):
    return json.dumps(value, sort_keys=True, indent=2, default=_json_value,
                      allow_nan=False).encode("utf-8")


def file_sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


@contextmanager
def _locked(path):
    """Thread and advisory process locks; never remove an owned lock file."""
    import fcntl
    key = str(Path(path).resolve())
    with _LOCKS_GUARD:
        lock = _LOCKS.setdefault(key, threading.RLock())
    with lock, Path(path).open("a+b") as handle:
        fcntl.flock(handle.fileno(), fcntl.LOCK_EX)
        try:
            yield
        finally:
            fcntl.flock(handle.fileno(), fcntl.LOCK_UN)


def _publish(staged, target):
    """Copy scratch into an output-side pending publication, then rename.

    Scratch and saved-output filesystems may differ. Only the final publication
    copy lives beside the output; numerical/download staging stays in scratch.
    """
    descriptor, temporary = tempfile.mkstemp(prefix=f".{target.name}.publishing-",
                                              dir=target.parent)
    temporary = Path(temporary)
    try:
        with os.fdopen(descriptor, "wb") as destination, Path(staged).open("rb") as source:
            shutil.copyfileobj(source, destination)
            destination.flush()
            os.fsync(destination.fileno())
        os.replace(temporary, target)
        directory_fd = os.open(target.parent, os.O_RDONLY)
        try:
            os.fsync(directory_fd)
        finally:
            os.close(directory_fd)
    finally:
        temporary.unlink(missing_ok=True)


def _atomic_json(target, value, scratch_dir):
    descriptor, temporary = tempfile.mkstemp(prefix="native-metadata-", suffix=".json",
                                              dir=scratch_dir)
    temporary = Path(temporary)
    try:
        with os.fdopen(descriptor, "wb") as handle:
            handle.write(_json_bytes(value))
            handle.flush()
            os.fsync(handle.fileno())
        _publish(temporary, target)
    finally:
        temporary.unlink(missing_ok=True)


@dataclass(frozen=True)
class NativeROMSSource:
    """A collection URL with an optional explicit actual time coordinate."""
    name: str
    dataset_url: str
    time_name: str | None = None
    engine: str | None = "pydap"
    collection: bool = True

    def url(self, valid_time):
        return self.dataset_url

    def identity(self):
        return {"name": self.name, "dataset_url": self.dataset_url,
                "time_name": self.time_name, "engine": self.engine,
                "collection": self.collection}


class UCSCNativeSource(NativeROMSSource):
    def __init__(self, url=UCSC_URL, *, engine="pydap"):
        super().__init__("ucsc_ccsra_native", url, "time", engine, True)


class WCOFSNativeSource(NativeROMSSource):
    def __init__(self, *, cycle=None, url_template=None, engine="pydap"):
        from .wcofs_download import WCOFSDownloader
        super().__init__("wcofs_native", url_template or WCOFSDownloader.roms_url_template,
                         "ocean_time", engine, False)
        object.__setattr__(self, "cycle", None if cycle is None else _utc(cycle))

    def url(self, valid_time):
        from .wcofs_download import WCOFSDownloader
        return WCOFSDownloader.source_url(valid_time, source="roms", cycle=self.cycle,
                                         url_template=self.dataset_url)

    def identity(self):
        return {**super().identity(), "cycle": None if self.cycle is None else str(self.cycle)}


def wcofs_mirror_url(nominal_url):
    """Map only an exact NCEI native nowcast to NOAA's identical public object.

    Forecasts, regulargrid products, altered date paths and other sources have no
    mapping. Never search for or substitute a different cycle/valid timestamp.
    """
    uri = urlsplit(str(nominal_url))
    if (uri.scheme != "https" or uri.netloc not in {"www.ncei.noaa.gov", "ncei.noaa.gov"}
            or uri.query or uri.fragment):
        return None
    path = re.fullmatch(r"/thredds/dodsC/model-wcofs-files/(\d{4})/(\d{2})/([^/]+)", uri.path)
    if not path:
        return None
    basename = path[3]
    legacy = re.fullmatch(r"nos\.wcofs\.fields\.n(003|006|009|012|015|018|021|024)\.(\d{8})\.t03z\.nc", basename)
    current = re.fullmatch(r"wcofs\.t03z\.(\d{8})\.fields\.n(003|006|009|012|015|018|021|024)\.nc", basename)
    if not legacy and not current:
        return None
    date = legacy[2] if legacy else current[1]
    try:
        pd.Timestamp(year=int(date[:4]), month=int(date[4:6]), day=int(date[6:]))
    except ValueError:
        return None
    if path[1] != date[:4] or path[2] != date[4:6]:
        return None
    return f"https://noaa-ofs-pds.s3.amazonaws.com/wcofs.{date}/{basename}"


class NativeHTTPMirrorOpener:
    """Verified identical-file mirror, staged locally before native subsetting.

    Download mode requires caller-approved ComputePathHandler scratch. The full
    source is removed after the raw Dataset closes; only the archive's selected
    region is retained. Byte-range mode is optional for diagnostics. Neither mode
    changes provider identity, source attributes, timestamps or numerical values.
    """
    def __init__(self, mirror_url_for, *, primary_engine="pydap", timeout=120,
                 primary_opener=None, mirror_opener=None, use_byte_range=False,
                 scratch_dir=None, http_session=None, local_opener=None):
        if not np.isfinite(timeout) or timeout <= 0:
            raise ValueError("Native HTTP timeout must be finite and positive.")
        self.timeout = timeout
        self.use_byte_range = bool(use_byte_range)
        self.scratch_dir = None if scratch_dir is None else Path(scratch_dir)
        if not self.use_byte_range and (self.scratch_dir is None or not self.scratch_dir.is_dir()):
            raise ValueError("Native HTTP staging requires an existing approved scratch_dir.")
        self.mirror_url_for = mirror_url_for
        self.primary_opener = primary_opener or (
            lambda url: open_raw(url, engine=primary_engine, timeout=timeout))
        self.mirror_opener = mirror_opener or (
            lambda url: open_raw(url, engine="netcdf4"))
        self.local_opener = local_opener or (
            lambda path: open_raw(path, engine="netcdf4"))
        self.http_session = http_session
        self._preferred_mirrors = set()
        self._state_lock = threading.Lock()

    @staticmethod
    def _https(url):
        uri = urlsplit(str(url))
        if uri.scheme != "https" or not uri.netloc or uri.fragment:
            raise ValueError("Native mirror retrieval must use verified HTTPS without a fragment.")
        return str(url)

    def on_failure(self, nominal_url, error):
        if not isinstance(error, MissingNativeTime) and self.mirror_url_for(nominal_url):
            with self._state_lock:
                self._preferred_mirrors.add(str(nominal_url))

    def _download_mirror(self, mirror_url):
        """Stream verified bytes into approved scratch, rejecting partial sources."""
        import requests
        session = self.http_session or requests.Session()
        owns_session = self.http_session is None
        verification = getattr(session, "verify", True)
        if verification is not True and not (isinstance(verification, (str, os.PathLike)) and str(verification)):
            if owns_session:
                session.close()
            raise ValueError("Native HTTP mirror certificate verification must remain enabled.")
        descriptor, temporary = tempfile.mkstemp(prefix="native-http-source-", suffix=".nc",
                                                  dir=self.scratch_dir)
        temporary = Path(temporary)
        try:
            current_url = self._https(mirror_url)
            with os.fdopen(descriptor, "wb") as destination:
                for _ in range(6):
                    # Follow redirects ourselves so an HTTPS-to-HTTP downgrade
                    # is rejected before any unverified request is dispatched.
                    with session.get(current_url, stream=True, timeout=self.timeout,
                                     verify=verification, allow_redirects=False,
                                     headers={"Accept-Encoding": "identity"}) as response:
                        retrieval_url = self._https(getattr(response, "url", None) or current_url)
                        if response.status_code in {301, 302, 303, 307, 308}:
                            location = response.headers.get("Location")
                            if not location:
                                raise OSError("Native HTTP redirect has no Location header.")
                            current_url = self._https(urljoin(retrieval_url, location))
                            continue
                        response.raise_for_status()
                        if response.status_code != 200:
                            raise OSError(f"Native HTTP full source requires status 200, got {response.status_code}.")
                        if response.headers.get("Content-Encoding", "identity").lower() not in {"", "identity"}:
                            raise OSError("Native HTTP full source must be transferred without content encoding.")
                        length = response.headers.get("Content-Length")
                        expected = None if length is None else int(length)
                        if expected is not None and expected < 0:
                            raise OSError("Native HTTP source has an invalid Content-Length.")
                        received = 0
                        for chunk in response.iter_content(chunk_size=1024 * 1024):
                            if chunk:
                                destination.write(chunk)
                                received += len(chunk)
                        if expected is not None and received != expected:
                            raise OSError(f"Incomplete native HTTP source: received {received}, expected {expected} bytes.")
                        destination.flush()
                        os.fsync(destination.fileno())
                        break
                else:
                    raise OSError("Native HTTP mirror exceeded five redirects.")
            with temporary.open("rb") as handle:
                signature = handle.read(8)
            if not (signature[:4] in {b"CDF\x01", b"CDF\x02", b"CDF\x05"}
                    or signature == b"\x89HDF\r\n\x1a\n"):
                raise OSError("Native HTTP source is not a NetCDF classic or HDF5 file.")
            return temporary, retrieval_url
        except BaseException:
            temporary.unlink(missing_ok=True)
            raise
        finally:
            if owns_session:
                session.close()

    @contextmanager
    def __call__(self, nominal_url):
        nominal_url = str(nominal_url)
        mirror_url = self.mirror_url_for(nominal_url)
        with self._state_lock:
            use_mirror = nominal_url in self._preferred_mirrors
        dataset, primary_error, staged_path = None, None, None
        if not use_mirror:
            try:
                dataset = self.primary_opener(nominal_url)
            except Exception as error:
                primary_error = error
                self.on_failure(nominal_url, error)
        retrieval_url = nominal_url
        try:
            if dataset is None:
                if not mirror_url:
                    if primary_error is not None:
                        raise primary_error
                    raise OSError(f"No exact native HTTP mirror for {nominal_url}.")
                mirror_url = self._https(mirror_url)
                if self.use_byte_range:
                    dataset = self.mirror_opener(mirror_url + "#mode=bytes")
                    retrieval_url = mirror_url
                else:
                    staged_path, retrieval_url = self._download_mirror(mirror_url)
                    dataset = self.local_opener(staged_path)
                dataset.encoding["native_archive_retrieval_url"] = retrieval_url
            else:
                dataset.encoding.setdefault("native_archive_retrieval_url", nominal_url)
            yield dataset
        except Exception as error:
            if retrieval_url == nominal_url:
                self.on_failure(nominal_url, error)
            raise
        finally:
            try:
                if dataset is not None:
                    dataset.close()
            finally:
                if staged_path is not None:
                    staged_path.unlink(missing_ok=True)


def make_wcofs_mirror_opener(*, primary_engine="pydap", timeout=120,
                              primary_opener=None, mirror_opener=None,
                              use_byte_range=False, scratch_dir=None,
                              http_session=None, local_opener=None):
    """Exact native NOAA S3 fallback for NativeROMSArchive(opener=...).

    Default full-file staging needs an existing caller-approved scratch_dir and
    verified HTTPS. The full source is deleted after selected arrays are loaded.
    use_byte_range=True offers NetCDF-C's slower experimental byte reader. Both
    preserve nominal request IDs and record actual retrieval provenance.
    """
    return NativeHTTPMirrorOpener(wcofs_mirror_url, primary_engine=primary_engine,
                                  timeout=timeout, primary_opener=primary_opener,
                                  mirror_opener=mirror_opener,
                                  use_byte_range=use_byte_range, scratch_dir=scratch_dir,
                                  http_session=http_session, local_opener=local_opener)


def open_raw(url, *, engine="pydap", timeout=120, session=None):
    options = {} if engine is None else {"engine": engine}
    if engine == "pydap":
        options["backend_kwargs"] = {"timeout": timeout}
        if session is not None:
            options["backend_kwargs"]["session"] = session
    return xr.open_dataset(url, decode_cf=False, decode_times=False,
                           mask_and_scale=False, decode_coords=False, cache=False,
                           create_default_indexes=False, **options)


def native_times(dataset, time_name=None):
    """Decode only the actual shared u/v time axis, never auxiliary ocean_time."""
    for name in ("u", "v"):
        if name not in dataset:
            raise ValueError(f"Missing native current variable {name}.")
    common = set(dataset.u.dims) & set(dataset.v.dims)
    if time_name is None:
        candidates = [name for name, variable in dataset.variables.items()
                      if variable.ndim == 1 and variable.dims[0] in common
                      and " since " in str(variable.attrs.get("units", ""))
                      and (name == variable.dims[0]
                           or variable.attrs.get("standard_name") == "time"
                           or variable.attrs.get("axis") == "T")]
        if len(candidates) != 1:
            raise ValueError(f"Ambiguous native current time axes {candidates}; specify time_name.")
        time_name = candidates[0]
    if time_name not in dataset or dataset[time_name].ndim != 1:
        raise ValueError("Native time must be a one-dimensional source variable.")
    dim = dataset[time_name].dims[0]
    if dim not in common:
        raise ValueError(f"{time_name} is auxiliary metadata, not the u/v scheduling axis.")
    variable = dataset[time_name]
    decoded = xr.decode_cf(xr.Dataset({time_name: variable.variable}), decode_times=True)
    values = decoded[time_name].values
    if (not np.issubdtype(values.dtype, np.datetime64) or not values.size
            or np.isnat(values).any()):
        raise ValueError("Native time must decode to finite Gregorian timestamps.")
    if np.any(np.diff(values) <= np.timedelta64(0, "ns")):
        raise ValueError("Native timestamps must be unique and strictly increasing.")
    return time_name, dim, values


def _bounds(bounds):
    result = np.asarray(bounds, dtype=float)
    if result.shape != (4,) or not np.isfinite(result).all():
        raise ValueError("bounds must be finite (west, east, south, north).")
    west, east, south, north = result
    if not (-180 <= west < east <= 180 and -90 <= south < north <= 90):
        raise ValueError("Use a non-wrapping WGS84 geographic bounding rectangle.")
    return tuple(map(float, result))


def native_window(dataset, bounds):
    """Choose nearest source nodes, retaining their contiguous C-grid envelope."""
    west, east, south, north = _bounds(bounds)
    if not {"lat_rho", "lon_rho"}.issubset(dataset.variables):
        raise ValueError("Native ROMS lat_rho/lon_rho are required.")
    lat = np.asarray(dataset.lat_rho.values, dtype=float)
    raw_lon = np.asarray(dataset.lon_rho.values, dtype=float)
    lon = (raw_lon + 180) % 360 - 180  # Selection only; stored coordinates are untouched.
    if (dataset.lat_rho.dims != ("eta_rho", "xi_rho")
            or dataset.lon_rho.dims != dataset.lat_rho.dims or lat.shape != lon.shape
            or lat.ndim != 2 or min(lat.shape) < 2):
        raise ValueError("RHO coordinates must share a two-dimensional eta_rho/xi_rho grid.")
    valid = np.isfinite(lat) & np.isfinite(lon)
    if not valid.any():
        raise ValueError("No finite native geographic nodes.")
    separable = (valid.all() and np.allclose(lat, lat[:, :1], atol=1e-9, rtol=0)
                 and np.allclose(lon, lon[:1, :], atol=1e-9, rtol=0))
    if separable:
        rows = [int(np.argmin(np.abs(lat[:, 0] - x))) for x in (south, north)]
        cols = [int(np.argmin(np.abs(lon[0] - x))) for x in (west, east)]
        method = "nearest_separable_boundary_nodes"
    else:
        from scipy.spatial import cKDTree
        indices = np.argwhere(valid)
        scale = np.cos(np.deg2rad((south + north) / 2))
        tree = cKDTree(np.column_stack((lon[valid] * scale, lat[valid])))
        count = 2 * max(lat.shape) + 1
        ys, xs = np.linspace(south, north, count), np.linspace(west, east, count)
        boundary = np.concatenate((np.column_stack((np.full(count, west), ys)),
                                   np.column_stack((np.full(count, east), ys)),
                                   np.column_stack((xs, np.full(count, south))),
                                   np.column_stack((xs, np.full(count, north)))))
        selected = indices[tree.query(boundary * [scale, 1])[1]]
        inside = np.argwhere(valid & (lon >= west) & (lon <= east)
                             & (lat >= south) & (lat <= north))
        selected = np.concatenate((selected, inside))
        rows, cols = selected[:, 0], selected[:, 1]
        method = "nearest_curvilinear_edge_nodes_and_interior_envelope"
    y0, y1, x0, x1 = int(min(rows)), int(max(rows)), int(min(cols)), int(max(cols))
    if y0 == y1 or x0 == x1:
        raise ValueError("The requested native window needs at least two RHO nodes per axis.")
    limits = {}
    for dim, size in dataset.sizes.items():
        if dim.startswith("eta_") or dim.startswith("xi_"):
            is_eta = dim.startswith("eta_")
            base = lat.shape[0 if is_eta else 1]
            low, high = (y0, y1) if is_eta else (x0, x1)
            if size not in (base, base - 1):
                raise ValueError(f"Unexpected native staggered dimension {dim}={size}.")
            stop = high + 1 if size == base else high
            limits[dim] = [low, stop - 1]
    indexers = {dim: slice(low, high + 1) for dim, (low, high) in limits.items()}
    selected_lat, selected_lon = lat[y0:y1 + 1, x0:x1 + 1], lon[y0:y1 + 1, x0:x1 + 1]
    actual = dict(west=float(np.nanmin(selected_lon)), east=float(np.nanmax(selected_lon)),
                  south=float(np.nanmin(selected_lat)), north=float(np.nanmax(selected_lat)))
    requested = dict(zip(("west", "east", "south", "north"), (west, east, south, north)))
    return indexers, {"method": method, "native_indices": limits,
                      "requested_extent": requested, "actual_extent": actual,
                      "boundary_offsets_degrees": {k: actual[k] - requested[k] for k in actual}}


def _dependencies(variable):
    attrs = variable.attrs
    names = set()
    for key in ("coordinates", "ancillary_variables", "bounds", "grid", "grid_mapping",
                "time", "_CoordinateAxes", "node_coordinates", "face_coordinates",
                "edge1_coordinates", "edge2_coordinates"):
        names.update(re.findall(r"[A-Za-z_][\w]*", str(attrs.get(key, ""))))
    for key in ("formula_terms", "cell_measures"):
        names.update(re.findall(r":\s*([A-Za-z_][\w]*)", str(attrs.get(key, ""))))
    return names


def _retained_variables(dataset, time_dim, time_name):
    absent = [name for name in REQUIRED_FIELDS if name not in dataset]
    if absent:
        raise ValueError(f"Missing native current interpretation fields: {absent}.")
    if any(any(key in dataset[name].encoding for key in ("scale_factor", "add_offset"))
           for name in CORE_FIELDS if name in dataset):
        raise ValueError("Source currents were already decoded; open with decode_cf=False.")
    names = {name for name, var in dataset.variables.items() if time_dim not in var.dims}
    names.update(name for name in CORE_FIELDS if name in dataset)
    names.update((time_name,))
    pending = list(names)
    while pending:
        name = pending.pop()
        references = _dependencies(dataset[name]) | set(dataset[name].dims)
        for reference in references:
            if reference in dataset.variables and reference not in names:
                names.add(reference)
                pending.append(reference)
            elif reference not in dataset.variables and reference not in dataset.dims:
                raise ValueError(f"Native field {name} references missing dependency {reference}.")
    return sorted(names)


def grid_sha256(dataset):
    digest = hashlib.sha256()
    for name in sorted(set(GRID_FIELDS) & set(dataset.variables)):
        var = dataset[name]
        array = np.asarray(var.values)
        if array.dtype.kind in "biufc":
            array = array.astype(array.dtype.newbyteorder("<"), copy=False)
        digest.update(_json_bytes({"name": name, "dims": list(var.dims),
                                   "dtype": str(array.dtype), "shape": list(array.shape)}))
        digest.update(array.tobytes(order="C"))
    return digest.hexdigest()


def subset_native(dataset, bounds, valid_time, *, time_name=None):
    """Select raw indices only: retain all sigma layers and encoded values."""
    name, time_dim, times = native_times(dataset, time_name)
    requested = _utc(valid_time).to_datetime64()
    indices = np.flatnonzero(times == requested)
    if len(indices) != 1:
        raise MissingNativeTime(f"Native source has no exact timestamp {_utc(valid_time)}.")
    spatial, metadata = native_window(dataset, bounds)
    spatial[time_dim] = slice(int(indices[0]), int(indices[0]) + 1)
    names = _retained_variables(dataset, time_dim, name)
    result = dataset[names].isel(spatial).copy()
    # Source backend/chunk storage hints cannot be reused for a smaller NetCDF.
    # Raw packing/fill attributes and dtypes stay on the variables themselves.
    for variable in result.variables.values():
        variable.encoding = {"_FillValue": None} if "_FillValue" not in variable.attrs else {}
    result.attrs = {k: v for k, v in result.attrs.items() if k != "_NCProperties"}
    metadata.update(time_name=name, time_dimension=time_dim,
                    source_time_index=int(indices[0]), actual_time=str(_utc(valid_time)),
                    variables=names, grid_sha256=grid_sha256(result))
    return result, metadata


def write_raw_netcdf(dataset, path):
    """Write raw types directly, bypassing CF packing/string transformations."""
    from netCDF4 import Dataset
    with Dataset(path, "w", format="NETCDF4") as output:
        for name, size in dataset.sizes.items():
            output.createDimension(name, size)
        output.setncatts(dataset.attrs)
        for name, variable in dataset.variables.items():
            array = np.asarray(variable.values)
            if array.dtype.kind not in "iufS":
                raise ValueError(f"Unsupported raw NetCDF type for {name}: {array.dtype}.")
            if array.dtype.kind == "S" and array.dtype.itemsize != 1:
                raise ValueError(f"Native CHAR metadata for {name} must retain single-character dtype.")
            attrs = dict(variable.attrs)
            fill = attrs.pop("_FillValue", False)
            target = output.createVariable(name, array.dtype, variable.dims, fill_value=fill)
            target.setncatts(attrs)
            target.set_auto_maskandscale(False)
            target[...] = array


def _same_raw(saved, expected):
    if dict(saved.sizes) != dict(expected.sizes) or set(saved.variables) != set(expected.variables):
        raise ValueError("Saved native dimensions/variables differ from the source subset.")
    for name in expected.variables:
        a, b = saved[name], expected[name]
        if a.dims != b.dims or a.dtype.newbyteorder("=") != b.dtype.newbyteorder("="):
            raise ValueError(f"Saved native dtype/dimensions changed for {name}.")
        try:
            np.testing.assert_array_equal(a.values, b.values)
            np.testing.assert_equal(a.attrs, b.attrs)
        except AssertionError as error:
            raise ValueError(f"Saved raw values/attributes changed for {name}.") from error
    try:
        np.testing.assert_equal(saved.attrs, expected.attrs)
    except AssertionError as error:
        raise ValueError("Saved global source metadata changed.") from error


class NativeROMSArchive:
    """One source-native regional archive, with immutable request and safe resume.

    out_dir and scratch_dir must already be approved handler-derived locations.
    An existing finalized file with conflicting provenance/checksum is preserved
    and raises ArchiveConflictError. A valid orphan from interrupted manifest
    publication is adopted without another source request.
    """
    def __init__(self, source, bounds, *, out_dir, scratch_dir, retries=3,
                 retry_seconds=1.0, opener=None):
        if isinstance(retries, bool) or not isinstance(retries, int) or retries < 1:
            raise ValueError("retries must be a positive integer.")
        if not np.isfinite(retry_seconds) or retry_seconds < 0:
            raise ValueError("retry_seconds must be finite and nonnegative.")
        self.source, self.bounds = source, _bounds(bounds)
        self.out_dir, self.scratch_dir = Path(out_dir).resolve(), Path(scratch_dir).resolve()
        if not self.scratch_dir.is_dir():
            raise FileNotFoundError(f"Approved ComputePathHandler scratch folder is missing: {self.scratch_dir}")
        self.out_dir.mkdir(parents=True, exist_ok=True)
        self.retries, self.retry_seconds = retries, float(retry_seconds)
        self.opener = opener or (lambda url: open_raw(url, engine=source.engine))
        self.request = {"format": "native-roms-v1", "source": source.identity(),
                        "bounds": list(self.bounds), "selection": "nearest_native_index_window",
                        "raw_values": True, "all_native_sigma_levels": True}
        self.request_id = hashlib.sha256(_json_bytes(self.request)).hexdigest()
        self.manifest_path = self.out_dir / "manifest.json"
        with _locked(self.out_dir / "archive.lock"):
            path = self.out_dir / "request.json"
            if path.exists():
                if json.loads(path.read_text()) != {**self.request, "request_id": self.request_id}:
                    raise ArchiveConflictError("Existing native archive has a different request identity.")
            else:
                _atomic_json(path, {**self.request, "request_id": self.request_id}, self.scratch_dir)
            self._manifest()

    def _manifest(self):
        if not self.manifest_path.exists():
            return {"request_id": self.request_id, "records": {}}
        manifest = json.loads(self.manifest_path.read_text())
        if manifest.get("request_id") != self.request_id or not isinstance(manifest.get("records"), dict):
            raise ArchiveConflictError("Native manifest request identity/schema differs.")
        return manifest

    def _record(self, stamp, record):
        with _locked(self.out_dir / "archive.lock"):
            manifest = self._manifest()
            manifest["records"][stamp] = record
            _atomic_json(self.manifest_path, manifest, self.scratch_dir)

    def inspect_times(self, start, end):
        """Return actual collection timestamps in [start,end), preserving gaps."""
        if not self.source.collection:
            raise ValueError("Per-time WCOFS files must be validated by save(); inspect_times requires a collection.")
        start, end = _utc(start), _utc(end)
        if start >= end:
            raise ValueError("The requested time interval must be increasing.")
        with NETCDF_LOCK, self.opener(self.source.url(start)) as dataset:
            _, _, values = native_times(dataset, self.source.time_name)
        selected = values[(values >= start.to_datetime64()) & (values < end.to_datetime64())]
        return pd.DatetimeIndex(selected)

    def _existing(self, target, stamp, url):
        with _locked(self.out_dir / "archive.lock"):
            previous = self._manifest()["records"].get(stamp)
        if not target.exists():
            return None
        sha = file_sha256(target)
        if previous and previous.get("status") == "saved" and previous.get("sha256") != sha:
            raise ArchiveConflictError(f"Saved native checksum differs; preserved file {target}.")
        try:
            with NETCDF_LOCK, open_raw(target, engine="netcdf4") as saved:
                if (saved.attrs.get("native_archive_request_id") != self.request_id
                        or saved.attrs.get("native_archive_source_url") != url):
                    raise ArchiveConflictError("Existing native file has different source/request metadata.")
                metadata = json.loads(saved.attrs["native_archive_subset"])
                _, _, values = native_times(saved, metadata["time_name"])
                if len(values) != 1 or values[0] != _utc(stamp).to_datetime64():
                    raise ArchiveConflictError("Existing native file has a different actual time.")
                if grid_sha256(saved) != metadata["grid_sha256"]:
                    raise ArchiveConflictError("Existing native grid hash differs.")
                if set(metadata["variables"]) != set(saved.variables):
                    raise ArchiveConflictError("Existing native interpretation fields differ.")
        except ArchiveConflictError:
            raise
        except Exception as error:
            raise ArchiveConflictError(f"Cannot validate existing native file {target}: {error}") from error
        record = {**metadata, "retrieval_url": metadata.get("retrieval_url", url),
                  "status": "saved", "url": url, "path": target.name,
                  "request_id": self.request_id, "sha256": sha,
                  "bytes": target.stat().st_size, "attempts": (previous or {}).get("attempts", 0),
                  "resumed": True, "error": ""}
        self._record(stamp, record)
        return record

    def save(self, valid_time):
        """Save exactly one native timestamp, or resume its verified final file."""
        valid_time = _utc(valid_time)
        stamp = str(valid_time)
        url = self.source.url(valid_time)
        label = valid_time.strftime("%Y%m%dT%H%M%S")
        fraction = valid_time.value % 1_000_000_000
        if fraction:
            label += f".{fraction:09d}"
        target = self.out_dir / f"native_{label}_{self.request_id[:12]}.nc"
        with _locked(target.with_suffix(".nc.lock")):
            existing = self._existing(target, stamp, url)
            if existing is not None:
                return existing
            errors = []
            for attempt in range(1, self.retries + 1):
                descriptor, temporary = tempfile.mkstemp(prefix="native-subset-", suffix=".nc",
                                                          dir=self.scratch_dir)
                os.close(descriptor)
                temporary = Path(temporary)
                try:
                    with NETCDF_LOCK, self.opener(url) as source:
                        subset, metadata = subset_native(source, self.bounds, valid_time,
                                                         time_name=self.source.time_name)
                        metadata["retrieval_url"] = source.encoding.get("native_archive_retrieval_url", url)
                        subset.attrs.update(native_archive_request_id=self.request_id,
                                            native_archive_source_url=url,
                                            native_archive_subset=json.dumps(metadata, sort_keys=True))
                        subset.load()
                    with NETCDF_LOCK:
                        write_raw_netcdf(subset, temporary)
                        with open_raw(temporary, engine="netcdf4") as saved:
                            _same_raw(saved, subset)
                    sha = file_sha256(temporary)
                    _publish(temporary, target)
                    record = {**metadata, "status": "saved", "url": url, "path": target.name,
                              "request_id": self.request_id, "sha256": sha,
                              "bytes": target.stat().st_size, "attempts": attempt,
                              "resumed": False, "error": "", "previous_errors": errors}
                    self._record(stamp, record)
                    return record
                except Exception as error:
                    errors.append(f"{type(error).__name__}: {error}")
                    failure_hook = getattr(self.opener, "on_failure", None)
                    if failure_hook is not None:
                        failure_hook(url, error)
                    if target.exists():
                        # Output publication succeeded but manifest publication failed.
                        # Preserve the orphan for verification/adoption on the next call.
                        raise
                    if isinstance(error, MissingNativeTime) or attempt == self.retries:
                        record = {"status": "missing" if isinstance(error, MissingNativeTime) else "failed",
                                  "url": url, "actual_time": None, "requested_time": stamp,
                                  "path": target.name, "request_id": self.request_id,
                                  "attempts": attempt, "error": errors[-1], "errors": errors}
                        self._record(stamp, record)
                        raise RuntimeError(f"Could not save native time {stamp}: {errors[-1]}") from error
                    if self.retry_seconds:
                        time.sleep(self.retry_seconds)
                finally:
                    temporary.unlink(missing_ok=True)
