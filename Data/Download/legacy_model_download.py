"""Shared chunk downloader for the legacy pydap velocity-model adapters.

The model object supplies its dataset, bounds, time indices, and filename
handler. Keeping those on the caller preserves existing regional archive paths
without depending on any plotting or HyperNav classes.
"""

import os
import pickle


def open_legacy_dataset(base_url, dataset_id, *, opener=None):
    """Open the configured pydap dataset, with an injectable legacy client."""
    if opener is None:
        from pydap.client import open_url
        opener = open_url
    return opener(base_url + dataset_id)


def open_legacy_copernicus_dataset(
    base_url, dataset_id, *, username, password,
    cas_url='https://cmems-cas.cls.fr/cas/login',
    opener=None, session_factory=None,
):
    """Keep legacy CAS session behavior separate from regional model classes."""
    if session_factory is None:
        from pydap.cas.get_cookies import setup_session
        session_factory = setup_session
    if opener is None:
        from pydap.client import open_url
        opener = open_url
    session = session_factory(cas_url, username, password)
    session.cookies.set('CASTGC', session.cookies.get_dict()['CASTGC'])
    return opener(base_url + dataset_id, session=session)


def _velocity_chunk(model, start, end, u_name, v_name, direct_data, time_as_list):
    selection = (
        slice(start, end),
        slice(None, len(model.depths)),
        slice(model.lllat_idx, model.urlat_idx),
        slice(model.lllon_idx, model.urlon_idx),
    )
    if direct_data:
        u = model.dataset[u_name].data[0][selection]
        v = model.dataset[v_name].data[0][selection]
        time = model.dataset['time'].data[start:end]
    else:
        u_grid = model.dataset[u_name][selection]
        v_grid = model.dataset[v_name][selection]
        u, v = u_grid[u_name].data, v_grid[v_name].data
        time = u_grid['time'].data
    if time_as_list:
        time = time.tolist()
    return {'u': u, 'v': v, 'time': time}


def download_legacy_chunks(
    model,
    *,
    u_name='water_u',
    v_name='water_v',
    recent=False,
    direct_data=False,
    time_as_list=False,
    reopen=None,
    print_index=False,
    create_parent=True,
):
    """Write the existing numbered pickle chunks, skipping cached files.

    ``recent`` retains the legacy ``len(indices) - 10`` starting index.
    ``reopen`` is called after a failed request and returns the replacement
    dataset. With no callback, errors propagate (the PACIOOS behavior).
    ``direct_data`` supports the older pydap array representation alongside
    grid variables; both write the same u/v/time archive dictionary.
    """
    indices = model.dataset_time.return_time_list()
    k = len(indices) - 10 if recent else 0
    while k < len(indices) - 1:
        filename = model.make_k_filename(k)
        print(k if print_index else filename)
        if create_parent:
            os.makedirs(os.path.dirname(filename), exist_ok=True)
        if os.path.isfile(filename):
            k += 1
            continue
        try:
            chunk = _velocity_chunk(
                model, indices[k], indices[k + 1],
                u_name, v_name, direct_data, time_as_list,
            )
            with open(filename, 'wb') as output:
                pickle.dump(chunk, output)
            k += 1
        except Exception:
            if reopen is None:
                raise
            print('Index ', k, ' encountered an error and did not save. Trying again')
            model.dataset = reopen()
