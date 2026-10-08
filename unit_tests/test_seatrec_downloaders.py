"""Offline regression tests for regional model transports and compatibility APIs."""

import importlib
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace
import unittest
from unittest.mock import Mock, patch

import numpy as np
import pandas as pd
import requests
import xarray as xr

from GeneralUtilities.Data.Download import copernicus_download as copernicus
from GeneralUtilities.Data.Download import hycom_forecast_download as forecast
from GeneralUtilities.Data.Download import hycom_october_download as archive
from GeneralUtilities.Data.Download._subset_download_utils import save_netcdf_atomic


def currents(times, *, depths=(0, 100, 600), variables=("water_u", "water_v"), value=1):
    times = pd.DatetimeIndex(times).tz_localize(None)
    shape = (len(times), len(depths), 2, 2)
    dataset = xr.Dataset(
        {name: (("time", "depth", "lat", "lon"), np.full(shape, value, dtype=float))
         for name in variables},
        coords={"time": times, "depth": list(depths), "lat": [17.0, 19.0], "lon": [-66.0, -64.0]},
    )
    dataset.depth.attrs.update(units="m", positive="down")
    for name in variables:
        dataset[name].attrs["units"] = "m/s"
    return dataset


def netcdf_bytes(dataset):
    with TemporaryDirectory(prefix="test_currents_") as folder:
        path = Path(folder) / "fixture.nc"
        dataset.to_netcdf(path, engine="netcdf4")
        return path.read_bytes()


class Response:
    def __init__(self, content, url="test://fixture"):
        self.content = content
        self.headers = {"Content-Length": str(len(content))}
        self.url = url
        self.ok = True
        self.status_code = 200

    def __enter__(self):
        return self

    def __exit__(self, *args):
        return False

    def raise_for_status(self):
        pass

    def iter_content(self, chunk_size):
        yield self.content


class ForecastSession:
    def __init__(self):
        self.headers = {}
        self.requests = []
        self.runs = pd.to_datetime(["2026-10-01T00:00Z", "2026-10-01T06:00Z"])

    def __enter__(self):
        return self

    def __exit__(self, *args):
        return False

    def get(self, url, params=None, **options):
        self.requests.append((url, params))
        if params is None:
            prefix = "FMRC_ESPC-D-V02_ts3z" if "ts3z" in url else "FMRC_ESPC-D-V02_uv3z"
            nodes = "".join(
                f'<dataset urlPath="{prefix}/runs/{prefix}_RUN_{run:%Y-%m-%dT%H:%M:%SZ}"/>'
                for run in self.runs
            )
            return Response(f"<catalog>{nodes}</catalog>".encode(), url)
        variables = params["var"].split(",")
        if "vertCoord" not in params:
            dataset = currents([self.runs[-1]], variables=variables)
        else:
            times = pd.date_range(params["time_start"], params["time_end"], freq="3h")
            dataset = currents(
                times, depths=[params["vertCoord"]], variables=variables,
                value=2 if "T06:00:00Z" in url else 1,
            )
        return Response(netcdf_bytes(dataset), url)


class ForecastTests(unittest.TestCase):
    def test_forecast_stitching_prefers_new_run_and_preserves_tracers(self):
        for tracer, expected in [(False, {"water_u", "water_v"}), (True, {"water_temp", "salinity"})]:
            with self.subTest(tracer=tracer):
                session = ForecastSession()
                with patch.object(forecast.requests, "Session", return_value=session):
                    dataset = forecast.NCSS_dataset(
                        -66, -64, 17, 19, "2026-10-01T03:00Z", "2026-10-01T09:00Z",
                        tracer=tracer, max_depth=100, max_attempts=1,
                    )
                self.assertEqual(set(dataset.data_vars), expected)
                np.testing.assert_array_equal(dataset.depth, [0, 100])
                np.testing.assert_array_equal(dataset.time, pd.date_range("2026-10-01T03:00", periods=3, freq="3h"))
                values = dataset[next(iter(expected))].isel(depth=0, lat=0, lon=0).values
                np.testing.assert_array_equal(values, [1, 2, 2])
                self.assertEqual(len(session.requests), 6)  # catalog, depth probe, four slabs

    def test_tds_retry_is_bounded_and_closes_failed_session(self):
        session = Mock()
        with patch.object(forecast.requests, "Session", return_value=session), patch.object(
            forecast.xr, "open_dataset", side_effect=OSError("transport unavailable")
        ) as opener:
            with self.assertRaisesRegex(OSError, "transport unavailable"):
                forecast.return_tds_dataset()
        self.assertEqual(opener.call_count, 2)
        session.close.assert_called_once()

    def test_compatibility_imports_share_forecast_functions(self):
        dataset = importlib.import_module("SeatrecPlots.Dataset")
        self.assertIs(dataset.NCSS_dataset, forecast.NCSS_dataset)
        self.assertIs(dataset.return_tds_dataset, forecast.return_tds_dataset)


class ArchiveTests(unittest.TestCase):
    def test_october_experiment_routing_includes_2017_transition(self):
        sources = archive.archive_sources(2017)
        self.assertEqual([url.rsplit("/", 2)[-2:] for url in sources["source_urls"]],
                         [["GLBv0.08", "expt_57.7"], ["expt_92.9", "uv3z"]])
        self.assertEqual(sources["variable_groups"], [["water_u", "water_v"]] * 2)
        self.assertIn("expt_53.X/data/2015", archive.archive_sources(2015)["source_urls"][0])
        self.assertIn("ESPC-D-V02/u3z/2024", archive.archive_sources(2024)["source_urls"][0])
        self.assertIs(archive.describe_hycom_source, archive.archive_sources)

    def test_archive_validation_preserves_masks_and_rejects_missing_times(self):
        times = pd.date_range("2024-10-01", periods=3, freq="3h", tz="UTC")
        dataset = currents(times)
        dataset.water_u.values[0, 0, 0, 0] = np.nan
        options = dict(expected_times=times, expected_depths=[0, 100, 600],
                       expected_grid={"lat": [17, 19], "lon": [-66, -64]},
                       variables=["water_u", "water_v"])
        result = archive._validate(dataset, **options)
        self.assertTrue(np.isnan(result.water_u.values[0, 0, 0, 0]))
        with self.assertRaisesRegex(ValueError, "timestamps"):
            archive._validate(dataset.isel(time=[0, 2]), **options)
        dataset.water_v.values[0, 0, 0, 0] = np.inf
        with self.assertRaisesRegex(ValueError, "infinity"):
            archive._validate(dataset, **options)

    def test_archive_publishes_complete_components_and_reuses_validated_output(self):
        times = pd.date_range("2024-10-01", periods=3, freq="3h", tz="UTC")
        sources = archive.archive_sources(2024)
        axes = {"time": times, "depth": np.array([0, 100, 600]),
                "lat": np.array([17, 19]), "lon": np.array([-66, -64]),
                "source_times": {url: times for url in sources["source_urls"]}}

        def subset(session, url, *, variables, times, depths, **kwargs):
            return currents(times, depths=depths, variables=variables)

        with TemporaryDirectory(prefix="test_archive_") as folder:
            options = dict(bounds={"west": -65.8, "east": -64.2, "south": 17.2, "north": 18.8},
                           max_depth=500, output_path=Path(folder) / "currents.nc",
                           cache_dir=Path(folder) / "cache", method="opendap", chunk_hours=6)
            with patch.object(archive, "_archive_axes", return_value=axes), patch.object(
                archive, "_opendap_subset", side_effect=subset
            ) as downloader:
                first = archive.download_hycom_october(times[0], times[-1], **options)
                second = archive.download_hycom_october(times[0], times[-1], **options)
            self.assertEqual(downloader.call_count, 2)  # one each for U and V
            np.testing.assert_array_equal(first.depth, [0, 100, 600])
            xr.testing.assert_identical(first, second)
            self.assertTrue(options["output_path"].is_file())
            self.assertIn("request_signature", first.attrs)

    def test_transport_failures_split_and_resume_cached_packets(self):
        times = pd.date_range("2024-10-01", periods=3, freq="3h", tz="UTC")

        def subset(session, url, *, variables, times, depths, **kwargs):
            if len(times) > 1:
                raise requests.Timeout("regional request timed out")
            return currents(times, depths=depths, variables=variables)

        with TemporaryDirectory(prefix="test_archive_packets_") as folder:
            options = dict(variables=["water_u"], times=times, depths=np.array([0, 100, 600]),
                           retained=np.array([0, 100]), grid={"lat": np.array([17, 19]), "lon": np.array([-66, -64])},
                           cache=Path(folder), signature_data={"processor": "test"}, signature="test-signature",
                           force=False, method="opendap", request_options={"timeout": 1})
            with patch.object(archive, "_opendap_subset", side_effect=subset) as downloader:
                first = archive._fetch_fragment(Mock(), "https://test/source", **options)
                calls = downloader.call_count
                second = archive._fetch_fragment(Mock(), "https://test/source", **options)
            self.assertEqual(calls, 5)  # two failures, then three successful single-time packets
            self.assertEqual(downloader.call_count, calls)
            xr.testing.assert_identical(first, second)

    def test_auto_transport_falls_back_to_ncss_and_retains_bracketing_depth(self):
        times = pd.date_range("2024-10-01", periods=3, freq="3h", tz="UTC")
        response = Response(netcdf_bytes(currents(times, variables=["water_u"])))
        session = Mock()
        session.get.return_value = response
        with TemporaryDirectory(prefix="test_archive_ncss_") as folder, patch.object(
            archive, "_opendap_subset", side_effect=ImportError("optional PyDAP absent")
        ):
            result = archive._fetch_fragment(
                session, "https://test/source", variables=["water_u"], times=times,
                depths=np.array([0, 100, 600]), retained=np.array([0, 100]),
                grid={"lat": np.array([17, 19]), "lon": np.array([-66, -64])},
                cache=Path(folder), signature_data={"processor": "test"},
                signature="test-signature", force=False, method="auto",
                request_options={"timeout": 1, "max_attempts": 1, "retry_wait": 0},
            )
        np.testing.assert_array_equal(result.depth, [0, 100])
        self.assertEqual(result.attrs["download_method"], "ncss")
        self.assertEqual(session.get.call_args.kwargs["timeout"], (20, 1))

    def test_truncated_ncss_response_retries_without_changing_retry_policy(self):
        truncated = Response(b"short")
        truncated.headers["Content-Length"] = "100"
        session = Mock()
        session.get.side_effect = [truncated, Response(b"complete")]
        with patch.object(archive.time, "sleep") as sleep:
            result = archive._request(session, "https://test/source", timeout=9,
                                      max_attempts=2, retry_wait=20, max_bytes=1000)
        self.assertEqual(result, b"complete")
        sleep.assert_called_once_with(15)
        self.assertEqual(session.get.call_count, 2)
        self.assertEqual(session.get.call_args.kwargs["timeout"], (20, 9))

    def test_october_compatibility_module_preserves_patchability(self):
        legacy = importlib.import_module("SeatrecPlots.october_currents")
        self.assertIs(legacy, archive)


class CopernicusTests(unittest.TestCase):
    def setUp(self):
        self.times = pd.date_range("2026-10-01", periods=3, freq="6h", tz="UTC")
        self.dataset = currents(self.times, variables=["uo", "vo"]).rename(
            {"lat": "latitude", "lon": "longitude"})
        self.bounds = {"west": -65.8, "east": -64.2, "south": 17.2, "north": 18.8}

    def sdk(self, data):
        def subset(**options):
            data.to_netcdf(Path(options["output_directory"]) / options["output_filename"], engine="netcdf4")
        return SimpleNamespace(subset=Mock(side_effect=subset))

    def test_sdk_options_loaded_result_and_consumer_validation_before_publish(self):
        sdk = self.sdk(self.dataset)
        with TemporaryDirectory(prefix="test_copernicus_") as folder:
            destination = Path(folder) / "currents.nc"
            destination.write_bytes(b"previous output")
            seen = []

            def validate(dataset):
                self.assertEqual(destination.read_bytes(), b"previous output")
                seen.append(dataset)

            with patch.object(copernicus.importlib, "import_module", return_value=sdk):
                result = copernicus.download_copernicus_data(
                    self.times[0], self.times[-1], bounds=self.bounds, output_path=destination,
                    username="user", password="secret", validate_dataset=validate,
                )
            self.assertEqual(len(seen), 1)
            xr.testing.assert_equal(result, self.dataset)
            options = sdk.subset.call_args.kwargs
            self.assertEqual(options["coordinates_selection_method"], "outside")
            self.assertEqual(options["variables"], ["uo", "vo"])
            self.assertEqual(options["password"], "secret")
            self.assertNotIn("password", result.attrs)
            with xr.open_dataset(destination, engine="netcdf4") as published:
                xr.testing.assert_equal(published, result)

    def test_invalid_coverage_and_consumer_failure_preserve_previous_output(self):
        for data, validator, message in [
            (self.dataset.isel(time=slice(1, None)), None, "time range"),
            (self.dataset.sel(depth=[0, 100]), None, "maximum depth"),
            (self.dataset, Mock(side_effect=ValueError("consumer rejected currents")), "consumer rejected"),
        ]:
            with self.subTest(message=message), TemporaryDirectory(prefix="test_copernicus_failure_") as folder:
                destination = Path(folder) / "currents.nc"
                destination.write_bytes(b"previous output")
                with patch.object(copernicus.importlib, "import_module", return_value=self.sdk(data)):
                    with self.assertRaisesRegex(ValueError, message):
                        copernicus.download_copernicus_data(
                            self.times[0], self.times[-1], bounds=self.bounds,
                            output_path=destination, validate_dataset=validator,
                        )
                self.assertEqual(destination.read_bytes(), b"previous output")

    def test_seatrec_adapter_keeps_waypoint_and_fieldset_contract(self):
        legacy = importlib.import_module("SeatrecPlots.copernicus_datasets")
        points = [{"lat": 18.0, "lon": -65.0}]
        fieldset = object()
        sdk = self.sdk(self.dataset)
        with patch.object(copernicus.importlib, "import_module", return_value=sdk), patch.object(
            legacy, "waypoint_bounds", return_value=self.bounds
        ), patch.object(legacy, "fieldset_from_dataset", return_value=fieldset):
            result, result_fieldset, result_points = legacy.download_copernicus_data(
                self.times[0], self.times[-1], waypoints=points,
            )
        xr.testing.assert_equal(result, self.dataset)
        self.assertIs(result_fieldset, fieldset)
        self.assertEqual(result_points, points)


class AtomicSaveTests(unittest.TestCase):
    def test_failed_write_preserves_output_and_removes_temporary_file(self):
        dataset = currents(pd.date_range("2024-10-01", periods=2, freq="3h"))
        with TemporaryDirectory(prefix="test_atomic_save_") as folder:
            path = Path(folder) / "currents.nc"
            path.write_bytes(b"previous output")
            with patch.object(xr.Dataset, "to_netcdf", side_effect=OSError("disk unavailable")):
                with self.assertRaisesRegex(OSError, "disk unavailable"):
                    save_netcdf_atomic(dataset, path, encode_time=True)
            self.assertEqual(path.read_bytes(), b"previous output")
            self.assertEqual(list(Path(folder).iterdir()), [path])


if __name__ == "__main__":
    unittest.main()
