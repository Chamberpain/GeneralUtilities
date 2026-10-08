"""Offline regressions for consolidated current downloaders and cache paths."""
import datetime
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from GeneralUtilities.Data.Download import hycom_download as hycom
from GeneralUtilities.Data.Download.download_paths import make_download_paths


def test_hycom_source_cache_is_shared_and_force_invalidates_metadata(monkeypatch):
    calls = []
    monkeypatch.setattr(hycom, 'open_url', lambda url, **kwargs: calls.append((url, kwargs)) or object())
    hycom.HYCOMDownloader._dataset_cache.clear()
    hycom.HYCOMDownloader._metadata_cache.clear()
    class Region(hycom.HYCOMDownloader):
        pass
    first = Region.get_dataset()
    assert Region.get_dataset() is first
    hycom.HYCOMDownloader._metadata_cache[Region.velocity_url] = {'old': True}
    assert Region.get_dataset(force=True) is not first
    assert len(calls) == 2
    assert Region.velocity_url not in Region._metadata_cache
    assert calls[0][1]['protocol'] == 'dap2'


def test_hycom_time_decode_preserves_legacy_class_and_utc():
    class LegacyTimeList(hycom.CurrentTimeList):
        pass
    class Region(hycom.HYCOMDownloader):
        time_list_class = LegacyTimeList
    variable = SimpleNamespace(attributes={'units': 'hours since 2024-01-01 00:00:00 UTC'})
    times, origin, step = Region._decode_time(variable, np.array([0, 3, 6]))
    assert isinstance(times, LegacyTimeList)
    assert step == 3600
    assert times[1] == origin + datetime.timedelta(hours=3)


def test_hycom_pickle_download_resume_and_corrupt_chunk_repair(tmp_path):
    class Region(hycom.HYCOMDownloader):
        location = 'test'
        file_handler = SimpleNamespace(tmp_file=lambda name: str(tmp_path / name))
        ref_date = datetime.datetime(2024, 1, 1)
        dataset_time = hycom.CurrentTimeList([datetime.datetime(2024,1,1) + datetime.timedelta(hours=i*3) for i in range(5)])
        depths = [-1, -10]
        lats, lons = [1, 2], [-2, -1]
        lllat_idx = lllon_idx = 0
        urlat_idx = urlon_idx = 2
        dataset = {'water_u': np.ones((5,2,2,2)), 'water_v': np.ones((5,2,2,2))*2}
        reads = 0
        @classmethod
        def initialize(cls, force=False):
            return cls
        @classmethod
        def _read_grid(cls, dataset, name, index):
            cls.reads += 1
            return dataset[name][index]
    result = Region.download_and_save(chunk_size=2, refresh_source=False)
    assert result['downloaded'] == [0,1,2]
    assert Region.reads == 6
    result = Region.download_and_save(chunk_size=2, refresh_source=False)
    assert result['skipped'] == [0,1,2]
    assert Region.reads == 6
    Path(Region.make_k_filename(1)).write_bytes(b'interrupted transfer')
    result = Region.download_and_save(chunk_size=2, refresh_source=False)
    assert result['redownloaded'] == [1]
    assert Region.reads == 8
    manifest = json.loads(Path(Region.manifest_filename()).read_text())
    assert all(row['status'] == 'complete' for row in manifest['chunks'].values())


def test_download_paths_create_only_tmp_and_out_and_reject_store(tmp_path, monkeypatch):
    import GeneralUtilities.Data.Filepath.instance as paths
    monkeypatch.setattr(paths, 'get_data_folder', lambda: str(tmp_path/'cache'))
    handler = make_download_paths(str(tmp_path/'Projects/Example/Utilities/Data'), 'HYCOM')
    assert Path(handler._tmp).is_dir()
    assert Path(handler._out).is_dir()
    assert not Path(handler._store).exists()
    assert not Path(handler._data).exists()
    with pytest.raises(PermissionError):
        handler.store_file('forbidden')


def test_hycom_dimensions_preserve_legacy_time_class(monkeypatch):
    class LegacyTimeList(hycom.CurrentTimeList):
        pass
    class Region(hycom.HYCOMDownloader):
        time_list_class = LegacyTimeList
    metadata = {'time': [datetime.datetime(2024,1,1)], 'lats': np.array([1.,2.]),
                'lons': np.array([-2.,-1.]), 'depths': np.array([-1.,-10.]),
                'units':'m/s', 'ref_date': datetime.datetime(2024,1,1)}
    monkeypatch.setattr(Region, '_source_metadata', classmethod(lambda cls,*args,**kwargs: metadata))
    dimensions = Region.get_dimensions(-1,-2,2,1,-10,{})
    assert isinstance(dimensions[0], LegacyTimeList)


def test_explicit_wcofs_builtin_processor_keeps_legacy_manifest_identity(monkeypatch):
    from GeneralUtilities.Data.Download import wcofs_download as wcofs
    from GeneralUtilities.Data.Download.wcofs_roms import regrid_roms
    class Region(wcofs.WCOFSDownloader):
        lllon, urlon, lllat, urlat, max_depth, location = -127,-123.9,38.9,43,700,'Humboldt'
    monkeypatch.setattr(wcofs,'download_archive',lambda fetch,**kwargs: kwargs['settings'])
    explicit = Region.download_record('2020-01-01','2020-01-02','unused',preprocess=regrid_roms)
    default = Region.download_record('2020-01-01','2020-01-02','unused')
    assert explicit['preprocess'] == default['preprocess'] == 'HyperNav.Utilities.Data.wcofs_roms:regrid_roms'
    assert explicit['preprocessor_sha256'] is None
    assert default['preprocessor_sha256'] is not None


def test_repair_loader_resolves_shared_directory_and_preserves_import_cache(monkeypatch, tmp_path):
    import sys
    from GeneralUtilities.Data.Download import hycom_redownload as redownload
    from GeneralUtilities.Data.Download import cec_hycom_download as repair
    sentinel = SimpleNamespace()
    monkeypatch.setitem(sys.modules,'hycom_download',sentinel)
    assert redownload.import_repair_function(Path(redownload.__file__).parent) is repair.repair_missing_files
    # The old project compatibility shim must also preserve an existing import.
    cec = tmp_path / 'CEC'
    cec.mkdir()
    (cec/'hycom_download.py').write_text(
        "import sys\n"
        "from GeneralUtilities.Data.Download import cec_hycom_download as _implementation\n"
        "sys.modules[__name__] = _implementation\n"
    )
    assert redownload.import_repair_function(cec) is repair.repair_missing_files
    assert sys.modules['hycom_download'] is sentinel
