"""Remote coordinates must be sliced before reads, including CF time decoding."""

import numpy as np
import pandas as pd
import pytest
import xarray as xr
from xarray.backends import BackendArray
from xarray.core import indexing

from GeneralUtilities.Data.Download import hycom_download as currents

NCSS = "https://ncss.hycom.org/thredds/ncss/grid/ESPC-D-V02/u3z/2024"
TIMES = pd.date_range("2024-10-01T03Z", "2024-10-01T09Z", freq="3h")
DEPTHS = np.array([0, 50, 100, 500])
GRID = {"lat": np.array([17.9, 18.0]), "lon": np.array([-65.0, -64.92])}


class LazyField(BackendArray):
    def __init__(self, name, values, reads):
        self.name = name
        self.values = values
        self.shape = values.shape
        self.dtype = values.dtype
        self.reads = reads

    def __getitem__(self, key):
        return indexing.explicit_indexing_adapter(
            key, self.shape, indexing.IndexingSupport.BASIC, self._read
        )

    def _read(self, key):
        self.reads.append((self.name, key))
        return self.values[key]


class Coordinate(LazyField):
    def _read(self, key):
        selected = np.arange(self.shape[0])[key[0]]
        if len(selected) == self.shape[0]:
            raise AssertionError(f"Unconstrained coordinate request: {self.name}")
        return super()._read(key)


def lazy_remote(monkeypatch, *, hourly=False, mismatch=False):
    reads, calls, closed = [], [], []
    times = pd.date_range(
        "2024-10-01",
        periods=36 if hourly else 12,
        freq="h" if hourly else "3h",
        tz="UTC",
    )
    axes = {
        "time": times,
        "depth": np.array([0, 50, 100, 500, 600]),
        "lat": np.array([17.8, 17.9, 18.0, 18.1]),
        "lon": np.array([0.0, 90.0, 180.0, 294.92, 295.0, 295.08, 295.16]),
    }
    coords = {}
    for name, axis in axes.items():
        attrs = {}
        if name == "time":
            values = (axis - times[0]).total_seconds().to_numpy() / 3600
            attrs = {
                "units": "hours since 2024-10-01 00:00:00",
                "calendar": "proleptic_gregorian",
            }
        else:
            values = axis
        coords[name] = xr.Variable(
            (name,),
            indexing.LazilyIndexedArray(Coordinate(name, values, reads)),
            attrs=attrs,
        )
    shape = tuple(len(axes[name]) for name in ("time", "depth", "lat", "lon"))
    values = np.arange(np.prod(shape), dtype=float).reshape(shape)
    source = xr.Dataset(
        {
            "water_u": xr.Variable(
                ("time", "depth", "lat", "lon"),
                indexing.LazilyIndexedArray(LazyField("water_u", values, reads)),
                attrs={"units": "m/s"},
            )
        },
        coords=xr.Coordinates(coords, indexes={}),
    )
    source.set_close(lambda: closed.append(True))

    def opened(url, **kwargs):
        calls.append((url, kwargs))
        assert kwargs["create_default_indexes"] is False
        assert kwargs["decode_times"] is False
        assert url.startswith("dap2://tds.hycom.org/")
        return source

    monkeypatch.setattr(currents.xr, "open_dataset", opened)
    if mismatch:
        axes = axes | {"time": times[:-1]}
    return source, axes, reads, calls, closed, values


@pytest.mark.parametrize("hourly", [False, True])
def test_no_full_axes_or_global_velocity_reads_and_time_decodes_exactly(
    monkeypatch, hourly
):
    source, axes, reads, calls, closed, values = lazy_remote(monkeypatch, hourly=hourly)
    loaded = currents.opendap_subset(
        object(),
        NCSS,
        variables=["water_u"],
        times=TIMES,
        depths=DEPTHS,
        grid=GRID,
        timeout=30,
        source_axes=axes,
    )
    assert len(calls) == 1 and closed == [True]
    assert not source.indexes
    np.testing.assert_allclose(loaded.lon, [295.0, 295.08])
    assert pd.DatetimeIndex(loaded.time.values).equals(TIMES.tz_localize(None))
    field = [(name, key) for name, key in reads if name == "water_u"]
    assert len(field) == 1
    selection = field[0][1]
    actual = tuple(
        part.indices(size) for part, size in zip(selection, values.shape, strict=True)
    )
    assert actual == (
        ((3, 10, 3) if hourly else (1, 4, 1)),
        (0, 4, 1),
        (1, 3, 1),
        (4, 6, 1),
    )
    np.testing.assert_array_equal(loaded.water_u.values, values[selection])
    for name, key in reads:
        if name != "water_u":
            assert len(np.arange(source.sizes[name])[key[0]]) < source.sizes[name]


def test_metadata_shape_mismatch_stops_before_any_remote_data_read(monkeypatch):
    _, axes, reads, _, closed, _ = lazy_remote(monkeypatch, mismatch=True)
    with pytest.raises(ValueError, match="time shape differs"):
        currents.opendap_subset(
            object(),
            NCSS,
            variables=["water_u"],
            times=TIMES,
            depths=DEPTHS,
            grid=GRID,
            timeout=30,
            source_axes=axes,
        )
    assert reads == [] and closed == [True]


def test_packed_velocities_are_scaled_once_and_missing_values_stay_nan(
    monkeypatch, tmp_path
):
    real_open = xr.open_dataset
    time = pd.date_range("2024-10-01", periods=12, freq="3h")
    z, y, x = (
        [0, 50, 100, 500, 600],
        [17.8, 17.9, 18.0, 18.1],
        [294.92, 295.0, 295.08, 295.16, 295.24],
    )
    data = np.full((12, 5, 4, 5), 1.234, dtype=np.float32)
    data[1, 2, 1, 1] = np.nan
    source = xr.Dataset(
        {"water_u": (("time", "depth", "lat", "lon"), data, {"units": "m/s"})},
        coords={"time": time, "depth": z, "lat": y, "lon": x},
    )
    path = tmp_path / "packed.nc"
    source.to_netcdf(
        path,
        encoding={
            "water_u": {"dtype": "int16", "scale_factor": 0.001, "_FillValue": -30000}
        },
    )

    def opened(url, **kwargs):
        kwargs["engine"] = "netcdf4"
        kwargs.pop("backend_kwargs")
        kwargs.pop("drop_variables")
        return real_open(path, **kwargs)

    monkeypatch.setattr(currents.xr, "open_dataset", opened)
    axes = {
        "time": pd.DatetimeIndex(time),
        "depth": np.array(z),
        "lat": np.array(y),
        "lon": np.array(x),
    }
    result = currents.opendap_subset(
        object(),
        NCSS,
        variables=["water_u"],
        times=TIMES,
        depths=DEPTHS,
        grid=GRID,
        timeout=30,
        source_axes=axes,
    )
    expected = source.isel(
        time=slice(1, 4), depth=slice(0, 4), lat=slice(1, 3), lon=slice(1, 3)
    )
    np.testing.assert_allclose(
        result.water_u, expected.water_u, rtol=0, atol=1e-6, equal_nan=True
    )
    assert np.isnan(result.water_u.values[0, 2, 0, 0])
    assert result.water_u.attrs["units"] == "m/s"
