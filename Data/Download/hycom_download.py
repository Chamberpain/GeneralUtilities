"""Constrained HYCOM OPeNDAP reads using verified source metadata.

Disable automatic coordinate indexes and postpone CF time decoding until the
requested slab is loaded. This prevents repeated whole-year coordinate reads.
The caller supplies its shared HTTP session and validates the returned grid.
"""

from urllib.parse import urlsplit

import numpy as np
import pandas as pd
import xarray as xr


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
