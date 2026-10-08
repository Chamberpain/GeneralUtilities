"""Convert retained WCOFS ROMS currents to a physical rectilinear grid.

The transform uses the file's h, zeta, s_rho, Cs_r, hc and Vtransform, following
https://www.myroms.org/wiki/Vertical_S-coordinate . C-grid velocities are
averaged to interior RHO points and rotated by the supplied grid angle.
Linear vertical interpolation uses depth below mean sea level, not sigma or
depth below the instantaneous surface. It never extrapolates above the top
RHO level or below the bottom RHO level. The outer native RHO row/column is
omitted because its velocity would require horizontal boundary extrapolation.

Horizontal interpolation uses neighboring native grid triangles, including
land vertices in the triangulation so masks and holes cannot be bridged. A
target gets a velocity only when every contributing vertex has both components.
No W is inferred from ROMS omega or terrain-following transport variables.
The default first fixed depth is 5 m; explicit shallower levels may be missing.
"""

from numbers import Real

import numpy as np
from scipy.spatial import Delaunay, QhullError
import xarray as xr


DEFAULT_DEPTHS = np.asarray(
    [
        5,
        10,
        15,
        20,
        25,
        30,
        40,
        50,
        60,
        75,
        100,
        125,
        150,
        200,
        250,
        300,
        400,
        500,
        600,
        700,
        800,
    ],
    dtype=float,
)
METRE_UNITS = {"m", "meter", "meters", "metre", "metres"}
VELOCITY_UNITS = {
    "m/s",
    "m s-1",
    "m s^-1",
    "m s**-1",
    "meter second-1",
    "meters second-1",
    "meters/second",
    "metres/second",
}


def _scalar(dataset, name):
    source = dataset[name] if name in dataset else dataset.attrs.get(name)
    if source is None or np.asarray(source).size != 1:
        raise ValueError(f"ROMS requires a scalar {name} parameter.")
    value = np.asarray(source).item()
    if isinstance(value, bool) or not isinstance(value, Real) or not np.isfinite(value):
        raise ValueError(f"ROMS {name} must be a finite number.")
    return float(value)


def regrid_roms(dataset, *, bounds, max_depth, depth_levels=None, grid_spacing=0.05):
    """Return loaded U/V(time, depth, lat, lon), in m/s and positive-down m.

    ``bounds`` is west, east, south, north. ``grid_spacing`` is a maximum
    spacing in degrees; both requested endpoints appear on each target axis.
    Supply physical ``depth_levels`` to control vertical resolution explicitly.
    Native fields must contain geographic RHO coordinates, radians ``angle``,
    metre ``h``/``zeta``, ``mask_rho``, and actual vertical transform metadata.
    Model water/dry masks and missing profiles remain missing in the result.
    """
    if not isinstance(dataset, xr.Dataset):
        raise TypeError("ROMS preprocessing requires an xarray.Dataset.")
    if (
        isinstance(grid_spacing, bool)
        or not isinstance(grid_spacing, Real)
        or not np.isfinite(grid_spacing)
        or grid_spacing <= 0
    ):
        raise ValueError("grid_spacing must be a finite positive number of degrees.")
    if (
        isinstance(max_depth, bool)
        or not isinstance(max_depth, Real)
        or not np.isfinite(max_depth)
        or max_depth <= 0
    ):
        raise ValueError("max_depth must be finite positive physical meters.")
    limits = np.asarray(bounds, dtype=float)
    if limits.shape != (4,) or not np.isfinite(limits).all():
        raise ValueError("bounds must contain finite west, east, south, north.")
    west, east, south, north = limits
    if (
        not (-180 <= west < east <= 180 and -90 <= south < north <= 90)
        or east - west > 180
    ):
        raise ValueError("Bounds must be an increasing regional geographic rectangle.")
    if depth_levels is None:
        levels = np.unique(
            np.append(DEFAULT_DEPTHS[DEFAULT_DEPTHS < max_depth], max_depth)
        )
    else:
        levels = np.asarray(depth_levels, dtype=float)
    if (
        levels.ndim != 1
        or not levels.size
        or not np.isfinite(levels).all()
        or (levels < 0).any()
        or (levels > max_depth).any()
        or (np.diff(levels) <= 0).any()
    ):
        raise ValueError(
            "depth_levels must increase uniquely from zero to max_depth, in meters."
        )

    required = {
        "u",
        "v",
        "lon_rho",
        "lat_rho",
        "angle",
        "h",
        "zeta",
        "s_rho",
        "Cs_r",
        "mask_rho",
    }
    missing = sorted(required - set(dataset.variables))
    if missing:
        raise ValueError(
            f"Missing ROMS geometry/transform/current variables: {missing}."
        )
    dataset = xr.decode_cf(dataset)
    time_name = "ocean_time" if "ocean_time" in dataset else "time"
    if time_name not in dataset or dataset[time_name].ndim != 1:
        raise ValueError("ROMS requires one-dimensional decoded ocean_time or time.")
    times = dataset[time_name].values
    if (
        not np.issubdtype(times.dtype, np.datetime64)
        or not len(times)
        or np.isnat(times).any()
        or (np.diff(times) <= np.timedelta64(0, "ns")).any()
    ):
        raise ValueError(
            "ROMS time must contain increasing nonmissing Gregorian timestamps."
        )
    time_dim = dataset[time_name].dims[0]
    lat, lon = dataset.lat_rho, dataset.lon_rho
    if lat.ndim != 2 or lat.dims != lon.dims or min(lat.shape) < 4:
        raise ValueError(
            "ROMS lon_rho/lat_rho must share a two-dimensional RHO grid of at least 4x4."
        )
    ydim, xdim = lat.dims
    for name in ("h", "angle", "mask_rho"):
        if dataset[name].dims != lat.dims:
            raise ValueError(f"ROMS {name} must use the RHO horizontal grid.")
    for name in ("h", "zeta"):
        if str(dataset[name].attrs.get("units", "")).lower().strip() not in METRE_UNITS:
            raise ValueError(f"ROMS {name} must declare physical metre units.")
    if str(dataset.angle.attrs.get("units", "")).lower().strip() not in {
        "rad",
        "radian",
        "radians",
    }:
        raise ValueError("ROMS angle must declare radians between XI and true east.")
    for name in ("lat_rho", "lon_rho"):
        expected = "degrees_north" if name == "lat_rho" else "degrees_east"
        if str(dataset[name].attrs.get("units", "")).lower().strip() not in {
            expected,
            expected.replace("degrees", "degree"),
        }:
            raise ValueError(f"ROMS {name} must declare geographic {expected}.")
    if dataset.zeta.dims != (time_dim, ydim, xdim):
        raise ValueError("ROMS zeta must use time and the RHO grid.")
    sigma = dataset.s_rho.values
    stretching = dataset.Cs_r.values
    if dataset.s_rho.dims != ("s_rho",) or dataset.Cs_r.dims != ("s_rho",):
        raise ValueError("s_rho and Cs_r must share the one-dimensional s_rho axis.")
    for name, values in (("s_rho", sigma), ("Cs_r", stretching)):
        if (
            len(values) < 2
            or not np.isfinite(values).all()
            or (values <= -1).any()
            or (values >= 0).any()
            or (np.diff(values) <= 0).any()
        ):
            raise ValueError(
                f"ROMS {name} must increase strictly between -1 and 0 at RHO levels."
            )
    transform, hc = _scalar(dataset, "Vtransform"), _scalar(dataset, "hc")
    if transform not in (1, 2) or hc < 0:
        raise ValueError("ROMS supports only Vtransform 1 or 2 and nonnegative hc.")
    if (
        "hc" in dataset
        and "units" in dataset.hc.attrs
        and str(dataset.hc.attrs["units"]).lower().strip() not in METRE_UNITS
    ):
        raise ValueError("ROMS hc must use physical metre units when declared.")

    # Read horizontal metadata first, then load just a regional rectangle and halo.
    lat_values = np.asarray(lat.values, dtype=float)
    lon_values = (np.asarray(lon.values, dtype=float) + 180) % 360 - 180
    region = (
        np.isfinite(lat_values)
        & np.isfinite(lon_values)
        & (lat_values >= south - grid_spacing)
        & (lat_values <= north + grid_spacing)
        & (lon_values >= west - grid_spacing)
        & (lon_values <= east + grid_spacing)
    )
    indices = np.argwhere(region)
    if not len(indices):
        raise ValueError("Requested geographic rectangle has no ROMS RHO points.")
    y0, x0 = np.maximum(indices.min(axis=0) - 2, 0)
    y1, x1 = np.minimum(indices.max(axis=0) + 3, lat.shape)
    if y1 - y0 < 4 or x1 - x0 < 4:
        raise ValueError(
            "Requested ROMS rectangle lacks enough neighboring velocity cells."
        )
    rho_indexers = {ydim: slice(y0, y1), xdim: slice(x0, x1)}
    h = np.asarray(dataset.h.isel(rho_indexers).values, dtype=float)
    angle = np.asarray(dataset.angle.isel(rho_indexers).values, dtype=float)
    mask = np.asarray(dataset.mask_rho.isel(rho_indexers).values)
    if not np.isfinite(mask).all() or not np.isin(mask, [0, 1]).all():
        raise ValueError("mask_rho must explicitly identify water=1 and land=0.")
    wet = mask == 1
    if (
        not wet.any()
        or not np.isfinite(h[wet]).all()
        or (h[wet] <= 0).any()
        or not np.isfinite(angle[wet]).all()
    ):
        raise ValueError(
            "Wet ROMS cells require finite positive bathymetry and finite angle."
        )
    if transform == 1 and hc > h[wet].min():
        raise ValueError(
            "Vtransform=1 requires hc no larger than minimum wet bathymetry."
        )
    h = np.where(wet, h, 1.0)
    angle = np.where(wet, angle, 0.0)
    zeta = np.asarray(dataset.zeta.isel(rho_indexers).values, dtype=float)
    wet_time = np.broadcast_to(wet, zeta.shape).copy()
    if "wetdry_mask_rho" in dataset:
        dry_mask = dataset.wetdry_mask_rho
        if dry_mask.dims != (time_dim, ydim, xdim):
            raise ValueError("wetdry_mask_rho must use time and the RHO grid.")
        dry_values = np.asarray(dry_mask.isel(rho_indexers).values)
        if not np.isfinite(dry_values).all() or not np.isin(dry_values, [0, 1]).all():
            raise ValueError("wetdry_mask_rho must contain finite zero/one masks.")
        wet_time &= dry_values == 1
    if (
        not np.isfinite(zeta[wet_time]).all()
        or (h[None, :, :] + zeta)[wet_time].min(initial=1) <= 0
    ):
        raise ValueError(
            "Wet ROMS columns require finite zeta and positive total water depth."
        )
    zeta = np.where(wet_time, zeta, 0.0)

    native = {}
    for name, expected_y, expected_x in (
        ("u", lat.shape[0], lat.shape[1] - 1),
        ("v", lat.shape[0] - 1, lat.shape[1]),
    ):
        variable = dataset[name]
        if (
            variable.ndim != 4
            or time_dim not in variable.dims
            or "s_rho" not in variable.dims
        ):
            raise ValueError(
                f"ROMS {name} must use time, s_rho and two staggered horizontal axes."
            )
        horizontal = [dim for dim in variable.dims if dim not in (time_dim, "s_rho")]
        uy, ux = horizontal
        if variable.sizes[uy] != expected_y or variable.sizes[ux] != expected_x:
            raise ValueError(
                f"ROMS {name} does not match the native Arakawa C-grid dimensions."
            )
        if str(variable.attrs.get("units", "")).lower().strip() not in VELOCITY_UNITS:
            raise ValueError(f"ROMS {name} must declare metre-per-second units.")
        indexers = {
            uy: slice(y0, y1 if name == "u" else y1 - 1),
            ux: slice(x0, x1 - 1 if name == "u" else x1),
        }
        values = np.asarray(
            variable.isel(indexers).transpose(time_dim, "s_rho", uy, ux).values,
            dtype=float,
        )
        if np.isinf(values).any():
            raise ValueError(
                f"ROMS {name} may contain NaN masks but no infinite velocities."
            )
        face_wet = (
            wet_time[:, :, :-1] & wet_time[:, :, 1:]
            if name == "u"
            else wet_time[:, :-1, :] & wet_time[:, 1:, :]
        )
        face_mask_name = "mask_u" if name == "u" else "mask_v"
        if face_mask_name in dataset:
            face_mask = dataset[face_mask_name]
            if face_mask.dims != (uy, ux):
                raise ValueError(
                    f"{face_mask_name} must use the corresponding staggered grid."
                )
            raw_mask = np.asarray(face_mask.isel(indexers).values)
            if not np.isfinite(raw_mask).all() or not np.isin(raw_mask, [0, 1]).all():
                raise ValueError(
                    f"{face_mask_name} must contain finite zero/one masks."
                )
            face_wet &= raw_mask[None, :, :] == 1
        moving_face_mask = f"wetdry_mask_{name}"
        if moving_face_mask in dataset:
            variable_mask = dataset[moving_face_mask]
            if variable_mask.dims != (time_dim, uy, ux):
                raise ValueError(
                    f"{moving_face_mask} must use time and its staggered grid."
                )
            raw_mask = np.asarray(variable_mask.isel(indexers).values)
            if not np.isfinite(raw_mask).all() or not np.isin(raw_mask, [0, 1]).all():
                raise ValueError(
                    f"{moving_face_mask} must contain finite zero/one masks."
                )
            face_wet &= raw_mask == 1
        native[name] = np.where(face_wet[:, None, :, :], values, np.nan)
    u = (native["u"][:, :, 1:-1, :-1] + native["u"][:, :, 1:-1, 1:]) / 2
    v = (native["v"][:, :, :-1, 1:-1] + native["v"][:, :, 1:, 1:-1]) / 2
    interior = (slice(1, -1), slice(1, -1))
    angle = angle[interior]
    eastward = u * np.cos(angle) - v * np.sin(angle)
    northward = u * np.sin(angle) + v * np.cos(angle)
    h = h[interior]
    eta = zeta[:, None, 1:-1, 1:-1]
    s, c = sigma[None, :, None, None], stretching[None, :, None, None]
    if transform == 1:
        z0 = hc * s + (h - hc) * c
        z = z0 + eta * (1 + z0 / h)
    else:
        z0 = (hc * s + h * c) / (hc + h)
        z = eta + (eta + h) * z0
    physical = -z[:, ::-1, :, :]
    wet_time = wet_time[:, 1:-1, 1:-1]
    if not (np.diff(physical, axis=1) > 0)[
        np.broadcast_to(wet_time[:, None, :, :], np.diff(physical, axis=1).shape)
    ].all():
        raise ValueError(
            "Computed ROMS physical depths must increase strictly down each wet column."
        )
    physical = physical.reshape(len(times), len(sigma), -1)
    eastward = eastward[:, ::-1].reshape(physical.shape)
    northward = northward[:, ::-1].reshape(physical.shape)
    wet_profiles = wet_time.reshape(len(times), -1)

    latitude = np.linspace(
        south, north, int(np.ceil((north - south) / grid_spacing - 1e-10)) + 1
    )
    longitude = np.linspace(
        west, east, int(np.ceil((east - west) / grid_spacing - 1e-10)) + 1
    )
    source_lat = lat_values[y0:y1, x0:x1][interior]
    source_lon = lon_values[y0:y1, x0:x1][interior]
    if (
        not np.isfinite(source_lat).all()
        or not np.isfinite(source_lon).all()
        or (np.abs(source_lat) > 90).any()
    ):
        raise ValueError(
            "The regional native RHO grid must have finite geographic coordinates."
        )
    scale = np.cos(np.deg2rad((south + north) / 2))
    source_points = np.column_stack((source_lon.ravel() * scale, source_lat.ravel()))
    if len(np.unique(source_points, axis=0)) != len(source_points):
        raise ValueError("ROMS geographic RHO coordinates must be unique.")
    target_lon, target_lat = np.meshgrid(longitude, latitude)
    queries = np.column_stack((target_lon.ravel() * scale, target_lat.ravel()))
    try:
        triangulation = Delaunay(source_points)
    except QhullError as exc:
        raise ValueError(
            "ROMS geographic points cannot form a valid horizontal grid."
        ) from exc
    simplex = triangulation.find_simplex(queries)
    valid_simplex = simplex >= 0
    triangles = triangulation.simplices[np.maximum(simplex, 0)]
    # Delaunay must not join nonneighboring cells across a folded grid or gap.
    rows, columns = np.unravel_index(triangles, source_lat.shape)
    valid_simplex &= np.ptp(rows, axis=1) <= 1
    valid_simplex &= np.ptp(columns, axis=1) <= 1
    transform_matrix = triangulation.transform[np.maximum(simplex, 0)]
    first_weights = np.einsum(
        "nij,nj->ni", transform_matrix[:, :2], queries - transform_matrix[:, 2]
    )
    weights = np.column_stack((first_weights, 1 - first_weights.sum(axis=1)))
    valid_simplex &= (weights >= -1e-10).all(axis=1)
    weights = np.maximum(weights, 0)
    weights /= weights.sum(axis=1, keepdims=True)
    contributing = weights > 1e-10
    output_shape = (len(times), len(levels), len(latitude), len(longitude))
    output_u, output_v = (
        np.full(output_shape, np.nan, dtype=np.float32) for _ in range(2)
    )
    for level_index, level in enumerate(levels):
        lower = np.clip((physical <= level).sum(axis=1) - 1, 0, len(sigma) - 2)
        indices = lower[:, None, :]
        d0 = np.take_along_axis(physical, indices, axis=1)[:, 0]
        d1 = np.take_along_axis(physical, indices + 1, axis=1)[:, 0]
        fraction = (level - d0) / (d1 - d0)
        valid_depth = (
            wet_profiles & (level >= physical[:, 0]) & (level <= physical[:, -1])
        )
        for native_values, output in ((eastward, output_u), (northward, output_v)):
            v0 = np.take_along_axis(native_values, indices, axis=1)[:, 0]
            v1 = np.take_along_axis(native_values, indices + 1, axis=1)[:, 0]
            available = (
                valid_depth
                & (np.isfinite(v0) | (fraction >= 1 - 1e-10))
                & (np.isfinite(v1) | (fraction <= 1e-10))
            )
            interpolated = (1 - fraction) * np.nan_to_num(
                v0, nan=0
            ) + fraction * np.nan_to_num(v1, nan=0)
            interpolated = np.where(available, interpolated, np.nan)
            gathered = interpolated[:, triangles]
            valid = valid_simplex[None, :] & (
                np.isfinite(gathered) | ~contributing[None, :, :]
            ).all(axis=2)
            values = (np.nan_to_num(gathered, nan=0) * weights[None, :, :]).sum(axis=2)
            output[:, level_index] = np.where(valid, values, np.nan).reshape(
                len(times), len(latitude), len(longitude)
            )
    jointly_valid = np.isfinite(output_u) & np.isfinite(output_v)
    if not jointly_valid.any():
        raise ValueError(
            "Requested physical grid has no jointly valid ROMS velocities; adjust bounds or depth levels."
        )
    output_u[~jointly_valid] = np.nan
    output_v[~jointly_valid] = np.nan
    return xr.Dataset(
        {
            "U": (
                ("time", "depth", "lat", "lon"),
                output_u,
                {"units": "m s-1", "standard_name": "eastward_sea_water_velocity"},
            ),
            "V": (
                ("time", "depth", "lat", "lon"),
                output_v,
                {"units": "m s-1", "standard_name": "northward_sea_water_velocity"},
            ),
        },
        coords={
            "time": ("time", times, {"standard_name": "time", "axis": "T"}),
            "depth": (
                "depth",
                levels,
                {
                    "units": "m",
                    "positive": "down",
                    "standard_name": "depth",
                    "axis": "Z",
                },
            ),
            "lat": (
                "lat",
                latitude,
                {"units": "degrees_north", "standard_name": "latitude", "axis": "Y"},
            ),
            "lon": (
                "lon",
                longitude,
                {"units": "degrees_east", "standard_name": "longitude", "axis": "X"},
            ),
        },
        attrs={
            "Conventions": "CF-1.8",
            "roms_vtransform": int(transform),
            "processing": "Interior RHO C-grid destaggering, angle rotation, linear physical-depth interpolation, neighboring geographic triangles; no extrapolation or mask bridging",
            "vertical_reference": "Depth below mean sea level",
            "maximum_grid_spacing_degrees": float(grid_spacing),
        },
    )
