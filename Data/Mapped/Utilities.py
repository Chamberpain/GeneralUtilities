"""Lazy rendering helpers for unchanged native mapped arrays."""


def pcolormesh_native(axis, longitude, latitude, values, **options):
    """Draw supplied coordinates/values without interpolation or grid conversion."""
    return axis.pcolormesh(longitude, latitude, values, **options)


def quiver_native(axis, longitude, latitude, eastward, northward, **options):
    """Draw supplied vector components at their original coordinates."""
    return axis.quiver(longitude, latitude, eastward, northward, **options)


def regional_axis(figure, rows=1, columns=1, index=1, *, map_template=None, projection=None):
    """Create a geographic axis with an optional caller-supplied map template."""
    import cartopy.crs as ccrs
    axis = figure.add_subplot(rows, columns, index, projection=projection or ccrs.PlateCarree())
    if map_template is not None:
        axis = map_template(ax=axis).get_map()[2]
    # Older map templates set removed xlabels_top/ylabels_right names.
    # Use current Gridliner flags so labels do not overlap adjacent panels.
    for artist in axis.get_children():
        if hasattr(artist, "top_labels") and hasattr(artist, "right_labels"):
            artist.top_labels = False
            artist.right_labels = False
    return axis
