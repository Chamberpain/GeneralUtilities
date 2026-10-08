"""Read and safely repair HYCOM files using validated manifest requests.

The physical reader and manifest repair are shared by model consumers.
"""

import hashlib
import os
from pathlib import Path
import re
import shutil
import tempfile
import time
from urllib.parse import parse_qs, urlsplit
import warnings

import numpy as np
import pandas as pd
import xarray as xr


def _read_uv(path):
    """Read physical velocities and normalize only coordinate conventions."""
    with xr.open_dataset(path, decode_times=True, mask_and_scale=True) as nc:
        frame = (
            nc[["water_u", "water_v"]].transpose("time", "depth", "lat", "lon").load()
        )
    if any(frame.sizes.get(c, 0) == 0 for c in ("time", "depth", "lat", "lon")):
        raise ValueError("Empty HYCOM time/depth/latitude/longitude axis.")
    if not np.issubdtype(frame.time.dtype, np.datetime64):
        raise ValueError("HYCOM time must decode to Gregorian datetime64.")
    frame = frame.assign_coords(lon=(frame.lon + 180.0) % 360.0 - 180.0)
    frame = frame.sortby(["time", "depth", "lat", "lon"])
    for name in ("depth", "lat", "lon"):
        a = frame[name].values
        if a.ndim != 1 or not np.isfinite(a).all() or np.any(np.diff(a) <= 0):
            raise ValueError(f"Invalid HYCOM {name} coordinates.")
    if str(frame.depth.attrs.get("positive", "down")).lower() != "down":
        raise ValueError("Expected positive-down HYCOM depths.")
    if str(frame.depth.attrs.get("units", "m")).lower() not in (
        "m",
        "meter",
        "meters",
        "metre",
        "metres",
    ):
        raise ValueError("Expected HYCOM depths in meters.")
    if pd.DatetimeIndex(frame.time.values).hasnans:
        raise ValueError("Missing HYCOM time coordinates.")
    for name in ("water_u", "water_v"):
        units = (
            str(frame[name].attrs.get("units", "m s-1"))
            .lower()
            .replace(" ", "")
            .replace("**", "^")
        )
        if units not in ("m/s", "ms-1", "ms^-1"):
            raise ValueError(f"{name} is not declared in m/s: {units!r}.")
    return frame


def repair_missing_files(
    paths,
    manifest_path=None,
    *,
    data_root=None,
    template_path=None,
    retries=4,
    timeout=(15, 180),
    backoff_seconds=2,
):
    """Download missing/empty manifest-listed chunks; return a list of repairs.

    manifest_path: optional; otherwise found by walking the paths' ancestors.
    data_root: base for relative output_path values (your CEC/Data directory).
    template_path: optional existing file with the SAME spatial/depth subset and
        time-of-day chunk as a missing file. Otherwise selected by product/part.
    retries: attempts per candidate URL. Permanent HTTP 4xx errors, other than
        408/429, immediately move to the next candidate.

    The first manifest write keeps download_manifest.csv.bak. Repaired files
    and CSV updates are installed with os.replace only after validation. CSV
    date, output_path, and unrelated columns/rows remain unchanged. Failures are
    recorded as download_failed and raised; days are never silently skipped.
    """
    paths = list(dict.fromkeys(Path(p).expanduser().resolve() for p in paths))
    missing = [p for p in paths if not p.is_file() or p.stat().st_size == 0]
    if not missing:
        return []
    if not isinstance(retries, int) or retries < 1 or backoff_seconds < 0:
        raise ValueError("retries must be >= 1 and backoff_seconds >= 0.")

    if manifest_path is None:
        found = {
            a / "download_manifest.csv"
            for p in paths
            for a in p.parents
            if (a / "download_manifest.csv").is_file()
        }
        if len(found) != 1:
            raise FileNotFoundError(
                f"HYCOM manifest refers to a missing file: {missing[0]}. "
                "Pass manifest_path=... (could not locate one unique manifest)."
            )
        manifest_path = found.pop()
    manifest_path = Path(manifest_path).expanduser().resolve()
    manifest = pd.read_csv(manifest_path, dtype=str, keep_default_na=False)
    if not {"date", "output_path"}.issubset(manifest.columns):
        raise ValueError("Manifest requires date and output_path columns.")
    roots = (
        [Path(data_root).expanduser().resolve()]
        if data_root
        else [manifest_path.parent.parent, manifest_path.parent]
    )
    resolved = []
    for raw in manifest.output_path:
        p = Path(raw).expanduser()
        candidates = (
            [p.resolve()] if p.is_absolute() else [(r / p).resolve() for r in roots]
        )
        matches = list(
            dict.fromkeys(p for p in candidates if p in paths or p.is_file())
        )
        if len(matches) > 1:
            raise ValueError(f"Ambiguous output_path {raw!r}; specify data_root.")
        resolved.append(matches[0] if matches else candidates[0])
    dates = (
        pd.to_datetime(manifest.date, utc=True, errors="raise")
        .dt.tz_convert(None)
        .dt.normalize()
    )
    if dates.isna().any():
        raise ValueError("Manifest contains missing dates.")
    rows_by_path = {}
    product_keys = []
    for i, path in enumerate(resolved):
        rows_by_path.setdefault(path, []).append(i)
        # The date and trailing request hash vary between otherwise matching chunks.
        product_key = re.sub(r"\d{4}-\d{2}-\d{2}", "DATE", path.name)
        product_keys.append(re.sub(r"_[0-9a-fA-F]{8,64}(?=\.nc$)", "", product_key))

    def commit(i, **values):
        # Re-read so other non-target rows/columns are not overwritten by an
        # earlier in-memory copy. Run one manifest writer at a time.
        fresh = pd.read_csv(manifest_path, dtype=str, keep_default_na=False)
        match = (fresh.output_path == manifest.at[i, "output_path"]) & (
            fresh.date == manifest.at[i, "date"]
        )
        if match.sum() != 1:
            raise ValueError(
                "Target manifest record was changed or duplicated during repair."
            )
        for key, value in values.items():
            if key not in fresh:
                fresh[key] = ""
            fresh.loc[match, key] = str(value)
        backup = manifest_path.with_suffix(manifest_path.suffix + ".bak")
        if not backup.exists():
            shutil.copy2(manifest_path, backup)
        fd, temp = tempfile.mkstemp(
            prefix=manifest_path.name + ".", suffix=".tmp", dir=manifest_path.parent
        )
        try:
            with os.fdopen(fd, "w", encoding="utf-8", newline="") as out:
                fresh.to_csv(out, index=False)
                out.flush()
                os.fsync(out.fileno())
            os.replace(temp, manifest_path)
        finally:
            Path(temp).unlink(missing_ok=True)

    reports = []
    import requests

    with requests.Session() as session:
        session.headers.update({"User-Agent": "CEC-HYCOM-missing-file-repair/1.0"})
        for target in missing:
            rows = rows_by_path.get(target, [])
            if len(rows) != 1:
                raise ValueError(
                    f"Expected one manifest row for {target}; found {len(rows)}."
                )
            i = rows[0]
            day = dates.iloc[i]
            row = manifest.iloc[i]
            last_url = ""
            try:
                if target.exists() and not target.is_file():
                    raise ValueError(f"Destination is not a regular file: {target}")
                named_date = re.search(r"\d{4}-\d{2}-\d{2}", target.name)
                if named_date and pd.Timestamp(named_date.group()) != day:
                    raise ValueError("Filename date disagrees with the manifest date.")
                # Reuse a recorded complete NCSS request, never an OPeNDAP
                # catalogue/HTML URL masquerading as a NetCDF download.
                urls = []
                for key in (
                    "repair_url",
                    "request_url",
                    "download_url",
                    "url",
                    "source_url",
                ):
                    url = str(row.get(key, "")).strip()
                    parsed = urlsplit(url)
                    q = parse_qs(parsed.query)
                    if (
                        parsed.scheme in ("http", "https")
                        and "/ncss/" in parsed.path
                        and {"north", "south", "east", "west", "var"}.issubset(q)
                        and ({"time_start", "time_end"}.issubset(q) or "time" in q)
                    ):
                        # Replay only a bounded, absolute request for this date,
                        # not a stale forecast URL or an entire-archive query.
                        try:
                            if "time" in q:
                                left = right = pd.to_datetime(
                                    q["time"][0], utc=True
                                ).tz_convert(None)
                            else:
                                left = pd.to_datetime(
                                    q["time_start"][0], utc=True
                                ).tz_convert(None)
                                right = pd.to_datetime(
                                    q["time_end"][0], utc=True
                                ).tz_convert(None)
                            if (
                                not pd.isna(left)
                                and not pd.isna(right)
                                and left <= right
                                and left >= day - pd.Timedelta(days=1)
                                and right <= day + pd.Timedelta(days=2)
                                and right - left <= pd.Timedelta(days=2)
                            ):
                                urls.append(url)
                        except (ValueError, TypeError):
                            pass
                # Prefer the closest intact chunk of the same product and
                # part number. The trailing hash can contain the request date.
                key = product_keys[i]
                candidates = []
                if template_path is not None:
                    p = Path(template_path).expanduser().resolve()
                    indices = rows_by_path.get(p, [])
                    m = re.search(r"\d{4}-\d{2}-\d{2}", p.name)
                    td = (
                        dates.iloc[indices[0]]
                        if indices
                        else (pd.Timestamp(m.group()) if m else None)
                    )
                    if td is None:
                        raise ValueError(
                            "Cannot determine template date; list it in the manifest or use a dated filename."
                        )
                    candidates = [(p, td)]
                else:
                    for j, p in enumerate(resolved):
                        if (
                            p != target
                            and product_keys[j] == key
                            and p.is_file()
                            and p.stat().st_size
                        ):
                            candidates.append((p, dates.iloc[j]))
                    candidates.sort(key=lambda item: abs(item[1] - day))
                template = None
                expected = None
                template_used = ""
                for p, td in candidates:
                    try:
                        template = _read_uv(p)
                        expected = pd.DatetimeIndex(template.time.values) + (day - td)
                        if expected.has_duplicates or expected[-1] - expected[
                            0
                        ] > pd.Timedelta(days=2):
                            raise ValueError(
                                "Template must be a unique daily or subdaily time chunk."
                            )
                        template_used = str(p)
                        break
                    except (OSError, ValueError, KeyError) as exc:
                        template = None
                        expected = None
                        if template_path is not None:
                            raise ValueError(
                                f"Invalid explicit template: {p}: {exc}"
                            ) from exc
                        warnings.warn(
                            f"Ignoring unusable repair template {p}: {exc}",
                            stacklevel=2,
                        )
                known_product = "gofs31a_expt93_0_glby" in target.name.lower()
                if template is not None and known_product:
                    west, east = float(template.lon.min()), float(template.lon.max())
                    south, north = float(template.lat.min()), float(template.lat.max())
                    if east - west >= 180:
                        raise ValueError(
                            "Automatic repair is intended for non-dateline-crossing regional subsets."
                        )
                    params = dict(
                        var="water_u,water_v",
                        west=west - 1e-4,
                        east=east + 1e-4,
                        south=south - 1e-4,
                        north=north + 1e-4,
                        horizStride=1,
                        timeStride=1,
                        time_start=expected[0].strftime("%Y-%m-%dT%H:%M:%SZ"),
                        time_end=expected[-1].strftime("%Y-%m-%dT%H:%M:%SZ"),
                        accept="netcdf",
                        addLatLon="true",
                    )
                    # Omitting vertCoord requests every native depth. Exact
                    # template levels/timestamps are selected after download.
                    base = "https://ncss.hycom.org/thredds/ncss/grid/GLBy0.08/expt_93.0"
                    endpoints = (
                        [base]
                        if expected[0].year != expected[-1].year
                        else [f"{base}/uv3z/{expected[0].year}", base]
                    )
                    urls.extend(
                        requests.Request("GET", u, params=params).prepare().url
                        for u in endpoints
                    )
                if not urls:
                    raise ValueError(
                        "Cannot reconstruct this subset safely. Supply an intact template_path "
                        "from the same product/part, or store the complete NCSS request_url "
                        "(variables, bounds and time) in this manifest row."
                    )
                # Without a template, a recorded GOFS request specifies its
                # own daily window. Validate against the native 3-hour grid.
                if expected is None:
                    q = parse_qs(urlsplit(urls[0]).query)
                    if "GLBy0.08/expt_93.0" not in urls[0]:
                        raise ValueError(
                            "A matching template_path is required for this product."
                        )
                    if "time" in q:
                        stamp = pd.to_datetime(
                            q["time"][0], utc=True, errors="raise"
                        ).tz_convert(None)
                        if stamp != stamp.floor("3h"):
                            raise ValueError(
                                "Recorded request must select a native 3-hour timestamp."
                            )
                        expected = pd.DatetimeIndex([stamp])
                    else:
                        left = (
                            pd.to_datetime(q["time_start"][0], utc=True)
                            .tz_convert(None)
                            .ceil("3h")
                        )
                        right = (
                            pd.to_datetime(q["time_end"][0], utc=True)
                            .tz_convert(None)
                            .floor("3h")
                        )
                        stride = int(q.get("timeStride", ["1"])[0])
                        if stride < 1:
                            raise ValueError("Invalid recorded timeStride.")
                        expected = pd.date_range(left, right, freq=f"{3 * stride}h")
                if (
                    len(expected) == 0
                    or expected[0] < day - pd.Timedelta(days=1)
                    or expected[-1] > day + pd.Timedelta(days=2)
                ):
                    raise ValueError(
                        "Request/template timestamps do not describe this manifest date."
                    )
                target.parent.mkdir(parents=True, exist_ok=True)
                errors = []
                repaired = False
                for url in dict.fromkeys(urls):
                    last_url = url
                    for attempt in range(1, retries + 1):
                        print(
                            f"Downloading missing HYCOM file: {target.name} (attempt {attempt}/{retries})",
                            flush=True,
                        )
                        fd, temp = tempfile.mkstemp(
                            prefix=target.stem + ".",
                            suffix=".part.nc",
                            dir=target.parent,
                        )
                        os.close(fd)
                        permanent_error = False
                        retry_after = 0.0
                        try:
                            with session.get(
                                url, stream=True, timeout=timeout
                            ) as response:
                                permanent_error = (
                                    400 <= response.status_code < 500
                                    and response.status_code not in (408, 429)
                                )
                                try:
                                    retry_after = max(
                                        0.0,
                                        float(response.headers.get("Retry-After", 0)),
                                    )
                                except ValueError:
                                    pass
                                response.raise_for_status()
                                with open(temp, "wb") as out:
                                    for chunk in response.iter_content(
                                        chunk_size=1024 * 1024
                                    ):
                                        if chunk:
                                            out.write(chunk)
                                    out.flush()
                                    os.fsync(out.fileno())
                            frame = _read_uv(
                                temp
                            )  # HTML, truncation, bad units etc. fail here.
                            if pd.DatetimeIndex(frame.time.values).has_duplicates:
                                raise ValueError(
                                    "Downloaded HYCOM timestamps contain duplicates."
                                )
                            actual = pd.DatetimeIndex(frame.time.values)
                            missing_times = expected.difference(actual)
                            if len(missing_times):
                                raise ValueError(
                                    f"Downloaded subset lacks required timestamps: {missing_times.tolist()}"
                                )
                            frame = frame.sel(time=expected)
                            if template is not None:
                                for name in ("depth", "lat", "lon"):
                                    # Roundoff tolerance only: NOT scientific spatial interpolation.
                                    idx = pd.Index(frame[name].values).get_indexer(
                                        template[name].values,
                                        method="nearest",
                                        tolerance=1e-5,
                                    )
                                    if (idx < 0).any() or len(np.unique(idx)) != len(
                                        idx
                                    ):
                                        raise ValueError(
                                            f"Downloaded {name} grid does not match the existing subset."
                                        )
                                    frame = frame.isel({name: idx}).assign_coords(
                                        {name: template[name].copy()}
                                    )
                            else:
                                q = parse_qs(urlsplit(url).query)
                                west, east = float(q["west"][0]), float(q["east"][0])
                                west, east = (
                                    (west + 180.0) % 360.0 - 180.0,
                                    (east + 180.0) % 360.0 - 180.0,
                                )
                                if not (
                                    west <= east
                                    and float(frame.lon.min()) >= west - 0.081
                                    and float(frame.lon.max()) <= east + 0.081
                                    and float(frame.lat.min())
                                    >= float(q["south"][0]) - 0.041
                                    and float(frame.lat.max())
                                    <= float(q["north"][0]) + 0.041
                                ):
                                    raise ValueError(
                                        "Downloaded domain does not match the recorded request bounds."
                                    )
                            good = np.isfinite(frame.water_u.values) & np.isfinite(
                                frame.water_v.values
                            )
                            if not good.any():
                                raise ValueError(
                                    "Downloaded subset has no finite paired U/V velocities."
                                )
                            # Save decoded physical velocities, not Parcels angular velocities.
                            for var in frame.variables:
                                frame[var].encoding = {}
                            frame.attrs.update(
                                repair_source_url=url, repair_template=template_used
                            )
                            frame.to_netcdf(temp)
                            with open(temp, "rb") as handle:
                                digest = hashlib.sha256()
                                for chunk in iter(
                                    lambda: handle.read(1024 * 1024), b""
                                ):
                                    digest.update(chunk)
                            nbytes = Path(temp).stat().st_size
                            with open(temp, "rb") as handle:
                                os.fsync(handle.fileno())
                            os.replace(temp, target)
                            repaired = True
                            break
                        except (
                            requests.RequestException,
                            OSError,
                            ValueError,
                            KeyError,
                            IndexError,
                        ) as exc:
                            errors.append(f"{url}: {type(exc).__name__}: {exc}")
                            if not permanent_error and attempt < retries:
                                time.sleep(
                                    max(
                                        retry_after,
                                        backoff_seconds * 2 ** (attempt - 1),
                                    )
                                )
                        finally:
                            Path(temp).unlink(missing_ok=True)
                        if permanent_error:
                            break
                    if repaired:
                        break
                if not repaired:
                    raise RuntimeError(
                        "No valid HYCOM replacement could be downloaded.\n"
                        + "\n".join(errors)
                    )
            except Exception as exc:
                commit(
                    i,
                    status="download_failed",
                    repair_status="failed",
                    repair_url=last_url,
                    repair_attempted_at_utc=pd.Timestamp.now(tz="UTC").isoformat(),
                    repair_error=str(exc),
                )
                raise RuntimeError(
                    f"Unable to repair {target}. The manifest records the failure. {exc}"
                ) from exc
            # Never label an unvalidated response as downloaded. If this CSV
            # commit fails, the validated file remains available for inspection.
            values = dict(
                status="downloaded",
                repair_status="repaired",
                repair_url=last_url,
                repaired_at_utc=pd.Timestamp.now(tz="UTC").isoformat(),
                file_size_bytes=nbytes,
                sha256=digest.hexdigest(),
                repair_template=template_used,
                repair_error="",
            )
            commit(i, **values)
            reports.append(dict(output_path=str(target), **values))
            print(
                f"Repaired {target.name}: {nbytes:,} bytes; manifest updated.",
                flush=True,
            )
    return reports
