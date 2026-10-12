"""Native ROMS geometry, encoded-value preservation and publication safety."""
from concurrent.futures import ThreadPoolExecutor
import json
from pathlib import Path
import threading

import numpy as np
import pandas as pd
import pytest
import xarray as xr

from GeneralUtilities.Data.Download import native_roms_archive as native


def dataset(*, curvilinear=False, auxiliary=False):
    ny, nx, nt, nz = 5, 6, 3, 3
    lat, lon = np.meshgrid(np.arange(40., 45.), np.arange(-127., -121.), indexing="ij")
    if curvilinear:
        lat = lat + 0.22 * np.sin(np.arange(nx))[None, :]
        lon = lon + 0.25 * np.sin(np.arange(ny))[:, None]
    actual = "time" if auxiliary else "ocean_time"
    ds = xr.Dataset(coords={actual: (actual, np.array([0., 24., 72.]),
                                      {"units": "hours since 2023-01-01 00:00:00",
                                       "calendar": "proleptic_gregorian"}),
                            "s_rho": ("s_rho", [-.8, -.5, -.2]),
                            "s_w": ("s_w", [-1., -.7, -.3, 0.])},
                    attrs={"Conventions": "CF-1.6 SGRID-0.3", "history": "unaltered source"})
    ds.s_rho.attrs.update(standard_name="ocean_s_coordinate_g2",
                          formula_terms="s: s_rho C: Cs_r eta: zeta depth: h depth_c: hc")
    ds.s_w.attrs.update(standard_name="ocean_s_coordinate_g2",
                        formula_terms="s: s_w C: Cs_w eta: zeta depth: h depth_c: hc")
    if auxiliary:
        ds["ocean_time"] = ("ocean_time", [12.], {"units": "seconds since 1940-01-01"})
    ds["time_run"] = (actual, np.array([0., 0., 48.]),
                       {"units": "hours since 2023-01-01"})
    grids = {"rho": (lat, lon), "u": ((lat[:, :-1] + lat[:, 1:]) / 2,
                                        (lon[:, :-1] + lon[:, 1:]) / 2),
             "v": ((lat[:-1] + lat[1:]) / 2, (lon[:-1] + lon[1:]) / 2),
             "psi": (lat[:-1, :-1], lon[:-1, :-1])}
    for grid, (ys, xs) in grids.items():
        dims = (f"eta_{grid}", f"xi_{grid}")
        ds[f"lat_{grid}"] = (dims, ys, {"units": "degrees_north"})
        ds[f"lon_{grid}"] = (dims, xs, {"units": "degrees_east"})
        ds[f"mask_{grid}"] = (dims, np.ones_like(xs, dtype="int8"))
    ds.mask_rho.values[0, 0] = 0
    for name, value in (("h", 500.), ("angle", .13), ("pm", .01), ("pn", .02)):
        ds[name] = (("eta_rho", "xi_rho"), np.full((ny, nx), value))
    for name, value in (("hc", 20.), ("Vtransform", 2), ("Vstretching", 4),
                        ("theta_s", 5.), ("theta_b", .4), ("Tcline", 20.)):
        ds[name] = xr.DataArray(value)
    ds["Cs_r"], ds["Cs_w"] = ds.s_rho.copy(), ds.s_w.copy()
    # Avoid self/eta dependencies on stretching arrays in this synthetic source.
    ds.Cs_r.attrs = {}; ds.Cs_w.attrs = {}
    ds["grid"] = xr.DataArray(np.int32(1), attrs={"cf_role": "grid_topology",
                                                  "node_coordinates": "lon_psi lat_psi"})
    for name, grid, vertical in (("u", "u", "s_rho"), ("v", "v", "s_rho"),
                                 ("w", "rho", "s_w")):
        shape = (nt, ds.sizes[vertical], *grids[grid][0].shape)
        raw = np.arange(np.prod(shape), dtype="int16").reshape(shape)
        raw[0, 0, 1, 1] = -32768
        ds[name] = ((actual, vertical, f"eta_{grid}", f"xi_{grid}"), raw,
                    {"units": "meter second-1", "scale_factor": np.float32(.01),
                     "add_offset": np.float32(.1), "_FillValue": np.int16(-32768),
                     "missing_value": np.int16(-32768), "grid": "grid",
                     "coordinates": f"lat_{grid} lon_{grid} {vertical} {actual} time_run",
                     "time": "ocean_time"})
    for name, grid in (("zeta", "rho"), ("ubar", "u"), ("vbar", "v")):
        ds[name] = ((actual, f"eta_{grid}", f"xi_{grid}"),
                    np.zeros((nt, *grids[grid][0].shape), dtype="float32"))
    # A non-current time-varying field must not be fetched/saved needlessly.
    ds["temp"] = ((actual, "s_rho", "eta_rho", "xi_rho"),
                   np.zeros((nt, nz, ny, nx), dtype="float32"))
    return ds


BOUNDS = (-126., -123., 41., 43.)


@pytest.fixture
def archive_factory(tmp_path):
    scratch = tmp_path / "scratch"
    scratch.mkdir()
    def create(ds=None, *, out="out", source=None, opener=None, **kwargs):
        ds = dataset() if ds is None else ds
        source = source or native.NativeROMSSource("fixture", "memory://fixture", engine=None)
        opener = opener or (lambda url: ds.copy(deep=True))
        return native.NativeROMSArchive(source, BOUNDS, out_dir=tmp_path / out,
                                         scratch_dir=scratch, retry_seconds=0,
                                         opener=opener, **kwargs)
    return create


def test_raw_packing_and_all_staggered_sigma_fields_roundtrip(archive_factory):
    ds = dataset()
    archive = archive_factory(ds)
    record = archive.save("2023-01-01")
    assert record["status"] == "saved" and not record["resumed"]
    assert record["native_indices"]["eta_rho"] == [1, 3]
    assert record["native_indices"]["xi_rho"] == [1, 4]
    path = archive.out_dir / record["path"]
    with native.open_raw(path, engine="netcdf4") as saved:
        assert saved.sizes["eta_u"] == 3 and saved.sizes["xi_u"] == 3
        assert saved.sizes["eta_v"] == 2 and saved.sizes["xi_v"] == 4
        assert saved.sizes["eta_psi"] == 2 and saved.sizes["xi_psi"] == 3
        assert saved.sizes["s_rho"] == 3 and saved.sizes["s_w"] == 4
        assert saved.u.dtype == np.dtype("int16")
        assert saved.u.values[0, 0, 0, 0] == -32768
        assert saved.u.attrs["scale_factor"] == ds.u.attrs["scale_factor"]
        for field in ("Cs_r", "Cs_w", "h", "angle", "pm", "pn", "hc", "grid",
                      "mask_rho", "mask_u", "mask_v", "mask_psi", "Vtransform"):
            assert field in saved
        np.testing.assert_array_equal(saved.u.values, ds.u.values[0:1, :, 1:4, 1:4])
        np.testing.assert_array_equal(saved.lat_rho.values, ds.lat_rho.values[1:4, 1:5])
        assert "temp" not in saved
        assert saved.attrs["history"] == ds.attrs["history"]
    assert record["sha256"] == native.file_sha256(path)
    assert not list(archive.scratch_dir.iterdir())


def test_true_collection_time_retains_auxiliary_ocean_time_and_gaps(archive_factory):
    ds = dataset(auxiliary=True)
    source = native.UCSCNativeSource(url="memory://ucsc", engine=None)
    archive = archive_factory(ds, source=source)
    times = archive.inspect_times("2023-01-01", "2023-01-05")
    assert list(times) == list(pd.to_datetime(["2023-01-01", "2023-01-02", "2023-01-04"]))
    record = archive.save("2023-01-04")
    with native.open_raw(archive.out_dir / record["path"], engine="netcdf4") as saved:
        assert saved.time.values.tolist() == [72.]
        assert saved.ocean_time.values.tolist() == [12.]
        assert saved.time_run.values.tolist() == [48.]
        assert saved.u.attrs["time"] == "ocean_time"
    with pytest.raises(ValueError, match="auxiliary"):
        native.native_times(ds, "ocean_time")


def test_curvilinear_nearest_envelope_slices_without_changing_values():
    ds = dataset(curvilinear=True)
    selected, metadata = native.subset_native(ds, BOUNDS, "2023-01-02")
    assert metadata["method"].startswith("nearest_curvilinear")
    for name in ("lat_rho", "lon_rho", "u", "v", "w", "mask_u"):
        indexers = {dim: slice(*[limits[0], limits[1] + 1])
                    for dim, limits in metadata["native_indices"].items() if dim in ds[name].dims}
        if "ocean_time" in ds[name].dims:
            indexers["ocean_time"] = slice(1, 2)
        np.testing.assert_array_equal(selected[name].values, ds[name].isel(indexers).values)
    assert selected.sizes["xi_u"] == selected.sizes["xi_rho"] - 1
    assert selected.sizes["eta_v"] == selected.sizes["eta_rho"] - 1


def test_native_longitudes_are_not_rewritten():
    ds = dataset()
    for name in ("lon_rho", "lon_u", "lon_v", "lon_psi"):
        ds[name].values[:] += 360
    selected, metadata = native.subset_native(ds, BOUNDS, "2023-01-01")
    assert selected.lon_rho.min() > 180
    assert metadata["actual_extent"]["west"] == -126


@pytest.mark.parametrize("missing", ["Cs_r", "hc", "angle", "grid"])
def test_interpretation_dependencies_cannot_be_silently_dropped(missing):
    ds = dataset().drop_vars(missing)
    with pytest.raises(ValueError, match="Missing|dependency"):
        native.subset_native(ds, BOUNDS, "2023-01-01")


def test_missing_native_time_records_missing_without_substitution(archive_factory):
    archive = archive_factory()
    with pytest.raises(RuntimeError, match="no exact timestamp"):
        archive.save("2023-01-03")
    manifest = json.loads(archive.manifest_path.read_text())
    record = manifest["records"]["2023-01-03 00:00:00"]
    assert record["status"] == "missing" and record["actual_time"] is None
    assert record["attempts"] == 1
    assert not list(archive.out_dir.glob("*.nc"))


def test_verified_resume_does_not_reopen_source(archive_factory):
    archive = archive_factory()
    first = archive.save("2023-01-01")
    archive.opener = lambda url: pytest.fail("Resume must not reopen the remote source")
    resumed = archive.save("2023-01-01")
    assert resumed["resumed"] and resumed["sha256"] == first["sha256"]


def test_request_or_checksum_conflicts_preserve_existing_files(archive_factory):
    archive = archive_factory()
    record = archive.save("2023-01-01")
    with pytest.raises(native.ArchiveConflictError, match="identity"):
        archive_factory(source=native.NativeROMSSource("different", "memory://fixture"))
    path = archive.out_dir / record["path"]
    with path.open("ab") as handle:
        handle.write(b"tampered")
    original = path.read_bytes()
    archive.opener = lambda url: pytest.fail("Conflicting finalized files must be preserved")
    with pytest.raises(native.ArchiveConflictError, match="checksum"):
        archive.save("2023-01-01")
    assert path.read_bytes() == original


def test_interrupted_manifest_publication_adopts_verified_orphan(archive_factory, monkeypatch):
    archive = archive_factory()
    original = archive._record
    monkeypatch.setattr(archive, "_record", lambda *args: (_ for _ in ()).throw(OSError("interrupted manifest")))
    with pytest.raises(OSError, match="interrupted"):
        archive.save("2023-01-01")
    assert len(list(archive.out_dir.glob("*.nc"))) == 1
    monkeypatch.setattr(archive, "_record", original)
    archive.opener = lambda url: pytest.fail("Valid orphan adoption must not redownload")
    adopted = archive.save("2023-01-01")
    assert adopted["resumed"] and adopted["attempts"] == 0
    assert json.loads(archive.manifest_path.read_text())["records"]


def test_failed_publication_cleans_scratch_and_retries_are_bounded(archive_factory, monkeypatch):
    archive = archive_factory(retries=2)
    original = native._publish
    def interrupt(staged, target):
        if target.suffix == ".nc":
            raise OSError("interrupted data publication")
        return original(staged, target)
    monkeypatch.setattr(native, "_publish", interrupt)
    with pytest.raises(RuntimeError, match="interrupted data publication"):
        archive.save("2023-01-01")
    assert not list(archive.out_dir.glob("*.nc"))
    assert not list(archive.scratch_dir.iterdir())
    record = json.loads(archive.manifest_path.read_text())["records"]["2023-01-01 00:00:00"]
    assert record["attempts"] == 2 and len(record["errors"]) == 2


def test_thread_lock_serializes_same_target_and_manifest_merge(archive_factory):
    count, guard = [0], threading.Lock()
    ds = dataset()
    def opener(url):
        with guard:
            count[0] += 1
        return ds.copy(deep=True)
    archive = archive_factory(opener=opener)
    with ThreadPoolExecutor(max_workers=3) as pool:
        results = list(pool.map(archive.save, ["2023-01-01", "2023-01-01", "2023-01-02"]))
    assert count[0] == 2
    assert sum(r["resumed"] for r in results) == 1
    assert len(json.loads(archive.manifest_path.read_text())["records"]) == 2


def test_wcofs_provider_reuses_native_source_routing():
    from GeneralUtilities.Data.Download.wcofs_download import WCOFSDownloader
    source = native.WCOFSNativeSource()
    stamp = "2023-02-26 00:00:00"
    assert source.url(stamp) == WCOFSDownloader.source_url(stamp, source="roms")
    assert ".fields." in source.url(stamp)


def test_decoded_packing_is_rejected_before_native_save():
    ds = dataset()
    ds.u.encoding["scale_factor"] = .01
    with pytest.raises(ValueError, match="already decoded"):
        native.subset_native(ds, BOUNDS, "2023-01-01")


def test_fractional_native_times_have_distinct_final_paths(archive_factory):
    ds = dataset()
    ds = ds.assign_coords(ocean_time=("ocean_time", [0., 1 / 3600, 1.5 / 3600],
                                     dict(ds.ocean_time.attrs)))
    archive = archive_factory(ds)
    first = archive.save("2023-01-01 00:00:01")
    second = archive.save("2023-01-01 00:00:01.5")
    assert first["path"] != second["path"]
    assert len(list(archive.out_dir.glob("*.nc"))) == 2


def test_native_char_scalars_and_arrays_keep_exact_dimensions(archive_factory):
    ds = dataset()
    ds["spherical"] = xr.DataArray(np.asarray(b"T", dtype="S1"))
    ds["native_label"] = ("label_char", np.asarray([b"a", b"b"], dtype="S1"))
    archive = archive_factory(ds)
    record = archive.save("2023-01-01")
    with native.open_raw(archive.out_dir / record["path"], engine="netcdf4") as saved:
        assert saved.spherical.dims == () and saved.spherical.dtype == np.dtype("S1")
        assert saved.native_label.dims == ("label_char",)
        np.testing.assert_array_equal(saved.native_label.values, ds.native_label.values)
    assert archive.save("2023-01-01")["resumed"]


def test_native_netcdf_das_namespace_attributes_remain_typed(archive_factory):
    # netcdf4 exposes native DAS containers as legal qualified NetCDF names.
    # Preserve those names/types rather than dropping or JSON-encoding metadata.
    ds = dataset()
    ds.attrs['EXTRA_DIMENSION.Nstate'] = np.int32(7)
    ds.attrs['DODS_EXTRA.Unlimited_Dimension'] = 'ocean_time'
    ds.u.attrs['_ChunkSizes'] = np.array([1, 3, 5, 5], dtype='int32')
    archive = archive_factory(ds)
    record = archive.save('2023-01-01')
    with native.open_raw(archive.out_dir / record['path'], engine='netcdf4') as saved:
        assert saved.attrs['EXTRA_DIMENSION.Nstate'] == np.int32(7)
        assert isinstance(saved.attrs['EXTRA_DIMENSION.Nstate'], np.int32)
        assert saved.attrs['DODS_EXTRA.Unlimited_Dimension'] == 'ocean_time'
        np.testing.assert_array_equal(saved.u.attrs['_ChunkSizes'], ds.u.attrs['_ChunkSizes'])
    assert archive.save('2023-01-01')['resumed']


LEGACY_NCEI = ("https://www.ncei.noaa.gov/thredds/dodsC/model-wcofs-files/2023/02/"
               "nos.wcofs.fields.n003.20230227.t03z.nc")
LEGACY_MIRROR = ("https://noaa-ofs-pds.s3.amazonaws.com/wcofs.20230227/"
                 "nos.wcofs.fields.n003.20230227.t03z.nc")


class HTTPResponse:
    def __init__(self, body, *, status=200, headers=None, url=None):
        self.body, self.status_code = body, status
        self.headers = {"Content-Length": str(len(body)), **(headers or {})}
        self.url = url
        self.closed = False

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.closed = True

    def raise_for_status(self):
        if self.status_code >= 400:
            raise OSError(f"HTTP {self.status_code}")

    def iter_content(self, chunk_size):
        # Include empty transport chunks and more than one populated chunk.
        yield self.body[:100]
        yield b""
        yield self.body[100:]


class HTTPSession:
    verify = True

    def __init__(self, *responses):
        self.responses = list(responses)
        self.calls = []

    def get(self, url, **kwargs):
        self.calls.append((url, kwargs))
        response = self.responses.pop(0)
        if response.url is None:
            response.url = url
        return response


def unavailable(url):
    raise OSError("primary archive unavailable")


def test_exact_mirror_mapping_preserves_native_basename_and_cycle():
    assert native.wcofs_mirror_url(LEGACY_NCEI) == LEGACY_MIRROR
    current = LEGACY_NCEI.replace("nos.wcofs.fields.n003.20230227.t03z.nc",
                                 "wcofs.t03z.20230227.fields.n003.nc")
    assert native.wcofs_mirror_url(current).endswith("wcofs.t03z.20230227.fields.n003.nc")
    for invalid in (LEGACY_NCEI.replace("n003", "f003"),
                    LEGACY_NCEI.replace(".fields.", ".fields.2ds."),
                    LEGACY_NCEI.replace("/2023/02/", "/2023/03/"),
                    LEGACY_NCEI.replace("20230227", "20230230"),
                    LEGACY_NCEI.replace("www.ncei.noaa.gov", "example.org"),
                    LEGACY_NCEI.replace("https:", "http:"), LEGACY_NCEI + "?alternate=true"):
        assert native.wcofs_mirror_url(invalid) is None


def test_primary_success_needs_no_mirror_and_preserves_source_attributes(tmp_path):
    scratch = tmp_path / "scratch"
    scratch.mkdir()
    source = dataset()
    session = HTTPSession()
    opener = native.make_wcofs_mirror_opener(scratch_dir=scratch,
                primary_opener=lambda url: source.copy(deep=True), http_session=session)
    with opener(LEGACY_NCEI) as opened:
        assert opened.attrs == source.attrs
        assert opened.encoding["native_archive_retrieval_url"] == LEGACY_NCEI
    assert not session.calls and not list(scratch.iterdir())


def test_full_mirror_staging_preserves_raw_subset_provenance_and_resume(archive_factory, tmp_path):
    source = native.WCOFSNativeSource(engine=None)
    ds = dataset()
    ds.ocean_time.attrs["units"] = "hours since 2023-02-26 06:00:00"
    payload = tmp_path / "fixture_source.nc"
    native.write_raw_netcdf(ds, payload)
    response = HTTPResponse(payload.read_bytes())
    session = HTTPSession(response)
    scratch = tmp_path / "scratch"
    events = []
    def local_opener(path):
        assert path.parent == scratch and path.exists()
        opened = native.open_raw(path, engine="netcdf4")
        original_close = opened._close
        def close():
            assert path.exists(), "Full source must remain until selected arrays load and source closes"
            events.append(("closed", path))
            original_close()
        opened.set_close(close)
        return opened
    # archive_factory creates the approved test scratch folder before opener setup.
    archive = archive_factory(ds, source=source)
    opener = native.make_wcofs_mirror_opener(primary_engine=None, scratch_dir=scratch,
                    primary_opener=unavailable, http_session=session, local_opener=local_opener)
    archive.opener = opener
    identity = archive.request_id
    record = archive.save("2023-02-26 06:00:00")
    assert record["url"] == LEGACY_NCEI and record["retrieval_url"] == LEGACY_MIRROR
    assert record["request_id"] == identity and record["attempts"] == 1
    assert session.calls[0][0] == LEGACY_MIRROR
    options = session.calls[0][1]
    assert options["verify"] is True and options["stream"] and not options["allow_redirects"]
    assert options["headers"]["Accept-Encoding"] == "identity"
    assert response.closed and len(events) == 1
    assert not events[0][1].exists() and not list(scratch.iterdir())
    with native.open_raw(archive.out_dir / record["path"], engine="netcdf4") as saved:
        assert saved.attrs["native_archive_source_url"] == LEGACY_NCEI
        metadata = json.loads(saved.attrs["native_archive_subset"])
        assert metadata["retrieval_url"] == LEGACY_MIRROR
        assert saved.u.dtype == ds.u.dtype
        np.testing.assert_array_equal(saved.u.values, ds.u.values[0:1, :, 1:4, 1:4])
        assert saved.attrs["history"] == ds.attrs["history"]
    archive.opener = lambda url: pytest.fail("Mirror resume must not fetch a full source again")
    resumed = archive.save("2023-02-26 06:00:00")
    assert resumed["resumed"] and resumed["retrieval_url"] == LEGACY_MIRROR
    assert resumed["request_id"] == identity


def test_byte_range_mode_is_optional_and_retains_bare_retrieval_url():
    calls = []
    def mirror(url):
        calls.append(url)
        return dataset()
    opener = native.make_wcofs_mirror_opener(use_byte_range=True,
                  primary_opener=unavailable, mirror_opener=mirror)
    with opener(LEGACY_NCEI) as opened:
        assert opened.encoding["native_archive_retrieval_url"] == LEGACY_MIRROR
    assert calls == [LEGACY_MIRROR + "#mode=bytes"]


@pytest.mark.parametrize("failure", ["truncated", "not_netcdf", "open", "read"])
def test_full_staging_failure_cleans_source_after_every_exit(tmp_path, failure):
    scratch = tmp_path / "scratch"
    scratch.mkdir()
    body = b"CDF\x02" + b"a" * 200
    headers = {"Content-Length": "999"} if failure == "truncated" else None
    if failure == "not_netcdf":
        body = b"<html>archive missing</html>"
    response = HTTPResponse(body, headers=headers)
    session = HTTPSession(response)
    paths = []
    def local(path):
        paths.append(path)
        if failure == "open":
            raise OSError("cannot open source")
        return dataset()
    opener = native.make_wcofs_mirror_opener(scratch_dir=scratch,
                  primary_opener=unavailable, http_session=session, local_opener=local)
    with pytest.raises(OSError):
        with opener(LEGACY_NCEI):
            if failure == "read":
                raise OSError("source payload read failed")
    assert response.closed and not list(scratch.iterdir())
    assert all(not path.exists() for path in paths)


def test_staging_rejects_disabled_tls_and_http_redirect_before_request(tmp_path):
    scratch = tmp_path / "scratch"
    scratch.mkdir()
    disabled = HTTPSession()
    disabled.verify = False
    opener = native.make_wcofs_mirror_opener(scratch_dir=scratch,
                  primary_opener=unavailable, http_session=disabled)
    with pytest.raises(ValueError, match="verification"):
        with opener(LEGACY_NCEI):
            pytest.fail("Unverified source must not open")
    assert not disabled.calls and not list(scratch.iterdir())
    redirect = HTTPResponse(b"", status=302, headers={"Location": "http://example.org/native.nc"})
    session = HTTPSession(redirect)
    opener = native.make_wcofs_mirror_opener(scratch_dir=scratch,
                  primary_opener=unavailable, http_session=session)
    with pytest.raises(ValueError, match="HTTPS"):
        with opener(LEGACY_NCEI):
            pytest.fail("HTTP downgrade must not open")
    assert len(session.calls) == 1 and redirect.closed and not list(scratch.iterdir())


def test_staging_requires_existing_approved_scratch(tmp_path):
    with pytest.raises(ValueError, match="approved scratch_dir"):
        native.make_wcofs_mirror_opener(scratch_dir=tmp_path / "unapproved")
    assert not (tmp_path / "unapproved").exists()


def test_primary_payload_failure_selects_exact_mirror_on_bounded_retry(archive_factory):
    ds = dataset()
    calls = []
    def primary(url):
        calls.append("primary")
        broken = ds.copy(deep=True)
        broken.attrs["nested_das"] = {"unserializable": 1}
        return broken
    def mirror(url):
        calls.append("mirror")
        return ds.copy(deep=True)
    opener = native.make_wcofs_mirror_opener(use_byte_range=True,
                   primary_opener=primary, mirror_opener=mirror)
    source = native.NativeROMSSource("fixture", LEGACY_NCEI, engine=None)
    archive = archive_factory(source=source, opener=opener, retries=2)
    record = archive.save("2023-01-01")
    assert record["attempts"] == 2 and calls == ["primary", "mirror"]
    assert record["retrieval_url"] == LEGACY_MIRROR
    assert not list(archive.scratch_dir.iterdir())


def test_mirror_does_not_substitute_wrong_native_timestamp(archive_factory, tmp_path):
    calls = []
    def mirror(url):
        calls.append(url)
        return dataset()
    opener = native.make_wcofs_mirror_opener(use_byte_range=True,
                  primary_opener=unavailable, mirror_opener=mirror)
    source = native.NativeROMSSource("fixture", LEGACY_NCEI, engine=None)
    archive = archive_factory(source=source, opener=opener, retries=3)
    with pytest.raises(RuntimeError, match="no exact timestamp"):
        archive.save("2023-01-03")
    assert len(calls) == 1 and not list(archive.out_dir.glob("*.nc"))
    assert not list(archive.scratch_dir.iterdir())


def test_resume_of_older_native_metadata_defaults_retrieval_to_primary(archive_factory):
    archive = archive_factory()
    record = archive.save("2023-01-01")
    target = archive.out_dir / record["path"]
    from netCDF4 import Dataset
    with Dataset(target, "a") as saved:
        metadata = json.loads(saved.native_archive_subset)
        del metadata["retrieval_url"]
        saved.native_archive_subset = json.dumps(metadata, sort_keys=True)
    manifest = json.loads(archive.manifest_path.read_text())
    previous = manifest["records"]["2023-01-01 00:00:00"]
    previous.pop("retrieval_url")
    previous["sha256"] = native.file_sha256(target)
    archive.manifest_path.write_text(json.dumps(manifest))
    archive.opener = lambda url: pytest.fail("Older request identity must remain resumable")
    resumed = archive.save("2023-01-01")
    assert resumed["resumed"] and resumed["retrieval_url"] == "memory://fixture"
