"""Offline archive checks without constructing model or /store directories."""

import pickle
from types import SimpleNamespace

import numpy as np
import pytest

from GeneralUtilities.Data.Download.legacy_model_download import (
    download_legacy_chunks, open_legacy_dataset, open_legacy_copernicus_dataset,
)
from GeneralUtilities.Data.Filepath.instance import FilePathHandler


class GridVariable:
    def __init__(self, name, values, times, *, failure=None):
        self.name, self.values, self.times = name, values, times
        self.failure = failure

    @property
    def data(self):
        if self.failure:
            raise self.failure
        return [self.values]

    def __getitem__(self, selection):
        if self.failure:
            raise self.failure
        return {
            self.name: SimpleNamespace(data=self.values[selection]),
            'time': SimpleNamespace(data=self.times[selection[0]]),
        }


def make_model(tmp_path, u_name='water_u', v_name='water_v', failure=None):
    values = np.arange(12 * 3 * 4 * 5).reshape(12, 3, 4, 5)
    times = np.arange(12) * 3
    dataset = {
        u_name: GridVariable(u_name, values, times, failure=failure),
        v_name: GridVariable(v_name, -values, times),
        'time': SimpleNamespace(data=times),
    }
    handler = object.__new__(FilePathHandler)
    handler._tmp = str(tmp_path)
    model = SimpleNamespace(
        dataset=dataset, depths=[0, -10],
        lllat_idx=1, urlat_idx=3, lllon_idx=2, urlon_idx=5,
        dataset_time=SimpleNamespace(return_time_list=lambda: list(range(12))),
        make_k_filename=lambda k: handler.tmp_file('GOFS_Hawaii_data/' + str(k)),
    )
    return model, values, times


@pytest.mark.parametrize('direct_data,time_as_list', [(False, False), (True, False), (True, True)])
def test_chunk_dictionary_and_existing_region_filename(tmp_path, direct_data, time_as_list):
    model, values, times = make_model(tmp_path)
    download_legacy_chunks(model, direct_data=direct_data, time_as_list=time_as_list)
    for k in range(11):
        with open(tmp_path / 'GOFS_Hawaii_data' / str(k), 'rb') as source:
            chunk = pickle.load(source)
        assert set(chunk) == {'u', 'v', 'time'}
        np.testing.assert_array_equal(chunk['u'], values[k:k + 1, :2, 1:3, 2:5])
        np.testing.assert_array_equal(chunk['v'], -values[k:k + 1, :2, 1:3, 2:5])
        np.testing.assert_array_equal(chunk['time'], times[k:k + 1])
        assert isinstance(chunk['time'], list) == time_as_list
    assert not (tmp_path / 'GOFS_Hawaii_data' / '11').exists()


def test_recent_download_preserves_chunk_numbers_and_skips_existing_files(tmp_path):
    model, _, _ = make_model(tmp_path)
    folder = tmp_path / 'GOFS_Hawaii_data'
    folder.mkdir()
    (folder / '3').write_bytes(b'existing archive')
    download_legacy_chunks(model, recent=True)
    assert sorted(int(path.name) for path in folder.iterdir()) == list(range(2, 11))
    assert (folder / '3').read_bytes() == b'existing archive'


def test_retry_reopens_and_replaces_dataset_before_retrying_same_chunk(tmp_path):
    model, _, _ = make_model(tmp_path, failure=OSError('unavailable'))
    healthy_model, values, _ = make_model(tmp_path)
    calls = []

    def reopen():
        calls.append('reopen')
        return healthy_model.dataset

    download_legacy_chunks(model, reopen=reopen)
    assert calls == ['reopen']
    assert model.dataset is healthy_model.dataset
    with open(tmp_path / 'GOFS_Hawaii_data' / '0', 'rb') as source:
        chunk = pickle.load(source)
    np.testing.assert_array_equal(chunk['u'], values[:1, :2, 1:3, 2:5])


def test_pacioos_grid_names_and_errors_without_reopen(tmp_path):
    model, values, _ = make_model(tmp_path, 'u', 'v')
    (tmp_path / 'GOFS_Hawaii_data').mkdir()
    download_legacy_chunks(model, u_name='u', v_name='v', print_index=True, create_parent=False)
    with open(tmp_path / 'GOFS_Hawaii_data' / '0', 'rb') as source:
        np.testing.assert_array_equal(pickle.load(source)['u'], values[:1, :2, 1:3, 2:5])

    broken, _, _ = make_model(tmp_path / 'broken', 'u', 'v', failure=OSError('unavailable'))
    with pytest.raises(OSError, match='unavailable'):
        download_legacy_chunks(broken, u_name='u', v_name='v')


def test_copernicus_direct_names_and_list_times(tmp_path):
    model, values, times = make_model(tmp_path, 'uo', 'vo')
    download_legacy_chunks(model, u_name='uo', v_name='vo', direct_data=True, time_as_list=True)
    with open(tmp_path / 'GOFS_Hawaii_data' / '0', 'rb') as source:
        chunk = pickle.load(source)
    np.testing.assert_array_equal(chunk['u'], values[:1, :2, 1:3, 2:5])
    assert chunk['time'] == times[:1].tolist()


def test_user_interrupt_is_not_retried(tmp_path):
    model, _, _ = make_model(tmp_path, failure=KeyboardInterrupt())
    with pytest.raises(KeyboardInterrupt):
        download_legacy_chunks(model, reopen=lambda: pytest.fail('must not retry interrupts'))


def test_legacy_provider_access_keeps_url_and_client_options():
    calls = []
    dataset = object()

    def opener(url, **options):
        calls.append((url, options))
        return dataset

    assert open_legacy_dataset('https://example.test/dodsC/', 'model', opener=opener) is dataset
    assert calls == [('https://example.test/dodsC/model', {})]


def test_copernicus_cas_cookie_and_session_are_forwarded():
    calls = []
    cookies = SimpleNamespace(
        get_dict=lambda: {'CASTGC': 'test-ticket'},
        set=lambda name, value: calls.append(('cookie', name, value)),
    )
    session = SimpleNamespace(cookies=cookies)
    dataset = object()

    def session_factory(url, username, password):
        calls.append(('session', url, username, password))
        return session

    def opener(url, **options):
        calls.append(('open', url, options))
        return dataset

    assert open_legacy_copernicus_dataset(
        'https://example.test/dodsC/', 'currents', username='test-user', password='test-password',
        opener=opener, session_factory=session_factory,
    ) is dataset
    assert calls == [
        ('session', 'https://cmems-cas.cls.fr/cas/login', 'test-user', 'test-password'),
        ('cookie', 'CASTGC', 'test-ticket'),
        ('open', 'https://example.test/dodsC/currents', {'session': session}),
    ]
