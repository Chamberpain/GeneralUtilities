"""Offline checks for forecast filenames, provider URLs, and pickle layouts."""

import importlib.util
import pickle
import sys
from datetime import datetime, timedelta
from pathlib import Path
from types import ModuleType, SimpleNamespace

import numpy as np
import pytest


@pytest.fixture
def downloader(monkeypatch):
    # Optional clients are replaced only within this fixture; no network access.
    pydap = ModuleType('pydap')
    client = ModuleType('pydap.client')
    client.open_url = lambda *args, **kwargs: pytest.fail('unexpected request')
    pydap.client = client
    webob = ModuleType('webob')
    exc = ModuleType('webob.exc')
    exc.HTTPError = type('HTTPError', (Exception,), {})
    webob.exc = exc
    for name, module in [('pydap', pydap), ('pydap.client', client), ('webob', webob), ('webob.exc', exc)]:
        monkeypatch.setitem(sys.modules, name, module)
    path = Path(__file__).parents[1] / 'Data/Download/production_hycom_download.py'
    name = 'GeneralUtilities.Data.Download.production_hycom_download'
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    monkeypatch.setitem(sys.modules, name, module)
    spec.loader.exec_module(module)

    class FixedDatetime(datetime):
        @classmethod
        def utcnow(cls):
            return cls(2026, 10, 8, 12)

    monkeypatch.setattr(module, 'datetime', FixedDatetime)
    return module


class VelocityVariable:
    def __init__(self, values):
        self.values = values
        self.reads = 0

    def __getitem__(self, selection):
        self.reads += 1
        return SimpleNamespace(data=self.values[selection])


class ForecastDataset(dict):
    attributes = {'NC_GLOBAL': {'time_origin': '2026-10-08 12:00:00'}}


def make_forecast():
    origin = datetime(2026, 10, 8, 12)
    times = np.array([(origin + timedelta(days=k)).timestamp() for k in range(8)])
    values = np.arange(8 * 8 * 9 * 11).reshape(8, 8, 9, 11)
    dataset = ForecastDataset(
        time=times, depth=np.arange(0, 800, 100),
        latitude=np.arange(15, 24), longitude=np.arange(200, 211),
        water_u=VelocityVariable(values), water_v=VelocityVariable(-values),
    )
    return dataset, values, times


def test_noaa_download_pickle_and_cached_forecast(downloader, monkeypatch, tmp_path):
    dataset, values, times = make_forecast()
    requests = []

    def open_url(url, **kwargs):
        requests.append((url, kwargs))
        return dataset

    monkeypatch.setattr(downloader, 'open_url', open_url)
    source = downloader.SourceNOAA()
    box = {'lllat': 16, 'urlat': 22, 'lllon': -159, 'urlon': -154}
    field, filename = source.download(box, str(tmp_path), 'HI')
    assert requests == [('https://www.ncei.noaa.gov/erddap/griddap/HYCOM_reg6_latest3d', {'output_grid': False})]
    assert Path(filename).name == 'HYCOM.NOAA.HI.2026100812Z.pickle'
    assert set(field) == {'data', 'dimensions'}
    np.testing.assert_array_equal(field['data']['U'], values[:5, :, 1:7, 1:6])
    np.testing.assert_array_equal(field['data']['V'], -values[:5, :, 1:7, 1:6])
    np.testing.assert_array_equal(field['data']['W'], np.zeros((5, 8, 6, 5)))
    np.testing.assert_array_equal(field['dimensions']['time'], times[:5])
    np.testing.assert_array_equal(field['dimensions']['lon'], np.arange(201, 206))
    with open(filename, 'rb') as output:
        stored = pickle.load(output)
    np.testing.assert_array_equal(stored['data']['U'], field['data']['U'])
    again, again_filename = source.download(box, str(tmp_path), 'HI')
    assert again_filename == filename
    np.testing.assert_array_equal(again['data']['U'], field['data']['U'])
    assert dataset['water_u'].reads == 1


def test_fmrc_run_fallback_and_origin(downloader, monkeypatch):
    calls = []

    def open_url(url, **kwargs):
        calls.append((url, kwargs))
        if len(calls) == 1:
            raise downloader.HTTPError('missing run')
        return {}

    monkeypatch.setattr(downloader, 'open_url', open_url)
    source = downloader.SourceFNMOCRun()
    assert source.get_datetime_origin() == datetime(2026, 10, 6, 12)
    assert calls[0][0].endswith('GLBy0.08_930_FMRC_RUN_2026-10-07T12:00:00Z')
    assert calls[1][0].endswith('GLBy0.08_930_FMRC_RUN_2026-10-06T12:00:00Z')
    assert all(kwargs == {'output_grid': False} for _, kwargs in calls)


def test_specific_fmrc_run_and_missing_run_error(downloader, monkeypatch):
    monkeypatch.setattr(downloader, 'open_url', lambda *args, **kwargs: {})
    source = downloader.SourceFNMOCRun('GLBy0.08_930_FMRC_RUN_2026-10-03T12:00:00Z')
    assert source.get_datetime_origin() == datetime(2026, 10, 3, 12)

    def missing(*args, **kwargs):
        raise downloader.HTTPError('missing run')

    monkeypatch.setattr(downloader, 'open_url', missing)
    with pytest.raises(downloader.HTTPError, match='requested'):
        downloader.SourceFNMOCRun('GLBy0.08_930_FMRC_RUN_2026-10-03T12:00:00Z')


def test_compatibility_module_is_same_shared_implementation(downloader, monkeypatch):
    import GeneralUtilities.Data.Download as package

    monkeypatch.setattr(package, 'production_hycom_download', downloader, raising=False)
    path = Path(__file__).parents[2] / 'HyperNav/Utilities/Production/hycom.py'
    if not path.exists():
        pytest.skip('legacy HyperNav checkout is not alongside GeneralUtilities')
    name = 'legacy_production_hycom'
    spec = importlib.util.spec_from_file_location(name, path)
    adapter = importlib.util.module_from_spec(spec)
    monkeypatch.setitem(sys.modules, name, adapter)
    spec.loader.exec_module(adapter)
    assert sys.modules[name] is downloader
