"""Native rendering must retain the supplied grid, values and masks."""
import numpy as np

from GeneralUtilities.Data.Mapped.Utilities import pcolormesh_native, quiver_native


class RecordingAxis:
    def pcolormesh(self, *args, **kwargs):
        return args, kwargs

    def quiver(self, *args, **kwargs):
        return args, kwargs


def test_pcolormesh_keeps_coordinate_and_value_identity():
    longitude = np.array([-127., -126.9])
    latitude = np.array([38.9, 39.])
    values = np.ma.masked_invalid([[np.nan, 1.], [2., 3.]])
    original = values.copy()
    args, options = pcolormesh_native(RecordingAxis(), longitude, latitude, values, shading="auto")
    assert args[0] is longitude and args[1] is latitude and args[2] is values
    assert options == {"shading": "auto"}
    np.testing.assert_array_equal(values.data, original.data)
    np.testing.assert_array_equal(values.mask, original.mask)


def test_quiver_keeps_native_components_without_rotation():
    longitude = np.array([-127., -126.9])
    latitude = np.array([38.9, 39.])
    eastward = np.array([[np.nan, 1.], [2., 3.]])
    northward = np.array([[4., 5.], [6., 7.]])
    originals = eastward.copy(), northward.copy()
    args, options = quiver_native(RecordingAxis(), longitude, latitude, eastward, northward, scale=7)
    assert args[0] is longitude and args[1] is latitude
    assert args[2] is eastward and args[3] is northward
    assert options == {"scale": 7}
    np.testing.assert_array_equal(eastward, originals[0])
    np.testing.assert_array_equal(northward, originals[1])
