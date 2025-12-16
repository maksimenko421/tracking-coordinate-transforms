import numpy as np
import random
import pytest

from geodesy import EarthEllipsoid


@pytest.fixture()
def get_latitudes():

    number_of_points = 100

    latitude2 = random.uniform(0, np.pi / 2)

    north_hemisphere = 1
    if north_hemisphere:
        latitude_start1 = 0
        latitude_end1 = np.pi / 2
    else:
        latitude_start1 = - np.pi / 2
        latitude_end1 = 0

    step_lat1 = (latitude_end1 - latitude_start1) / number_of_points

    latitudes1 = [latitude_start1 + step_lat1 * k for k in range(number_of_points)]

    return latitudes1, latitude2


def test_meridian_from_equator_to_meters_and_back(get_latitudes):
    # Checking meridian_to_meters_from_equator and meters_to_meridian_from_equator methods

    geo_eps = 1e-10
    meters_eps = 1e-3

    latitudes1, latitude2 = get_latitudes
    number_of_points = len(latitudes1)

    err_meters = np.zeros([number_of_points, ])
    err_latitude = np.zeros([number_of_points, ])

    for k in range(number_of_points):
        current_latitude = latitudes1[k]

        meridian_length = EarthEllipsoid.meridian_to_meters_from_equator(current_latitude)
        recalculate_lat_diff = EarthEllipsoid.meters_to_meridian_from_equator(meridian_length)
        if current_latitude < 0:
            recalculate_lat_diff = -recalculate_lat_diff
        recalculated_meridian_length = EarthEllipsoid.meridian_to_meters_from_equator(recalculate_lat_diff)

        err_meters[k] = meridian_length - recalculated_meridian_length
        err_latitude[k] = current_latitude - recalculate_lat_diff

    max_latitude_err = np.max(np.abs(np.array(err_latitude)))
    max_meters_err = np.max(np.abs(np.array(err_meters)))

    assert max_meters_err < meters_eps, "difference in meters is too big"
    assert max_latitude_err < geo_eps, "difference in radians is too big"


def test_meridian_to_meters_and_back(get_latitudes):
    # Checking meridian_to_meters and meters_to_meridian methods
 
    geo_eps = 1e-10
    meters_eps = 1e-3

    latitudes1, latitude2 = get_latitudes
    number_of_points = len(latitudes1)

    err_meters = np.zeros([number_of_points, ])
    err_latitude = np.zeros([number_of_points, ])
    for k in range(number_of_points):
        current_latitude1 = latitudes1[k]
        current_latitude2 = latitude2

        current_lat_diff = current_latitude2 - current_latitude1

        meridian_length = EarthEllipsoid.meridian_to_meters_with_sign(current_latitude1, current_latitude2)
        lat_diff = EarthEllipsoid.meters_to_meridian_with_sign(current_latitude1, meridian_length)

        recalculated_meridian_length = EarthEllipsoid.meridian_to_meters_with_sign(current_latitude1, current_latitude1 + lat_diff)

        err_meters[k] = np.abs(meridian_length - recalculated_meridian_length)
        err_latitude[k] = np.abs(current_lat_diff - lat_diff)

    max_latitude_err = np.max(np.abs(np.array(err_latitude)))
    max_meters_err = np.max(np.abs(np.array(err_meters)))

    assert max_meters_err < meters_eps, "difference in meters is too big"
    assert max_latitude_err < geo_eps, "difference in radians is too big"


def test_parallel_to_meters_and_back(get_latitudes):
    # Checking parallel_to_meters and meters_to_parallel methods

    geo_eps = 1e-10
    meters_eps = 1e-3

    sign = (2 * random.randint(0, 1) - 1)
    current_latitude = random.uniform(0, np.pi / 2) * sign

    number_of_points = 100

    lon_diff_start = 0.0
    lon_diff_end = 2 * np.pi
    step_lon = (lon_diff_end - lon_diff_start) / number_of_points

    lon_diffs = [lon_diff_start + step_lon * k for k in range(number_of_points)]

    err_meters = np.zeros([number_of_points, ])
    err_longitude = np.zeros([number_of_points, ])

    for k in range(number_of_points):

        current_lon_diff = lon_diffs[k]

        parallel_length = EarthEllipsoid.parallel_arc_length_with_sign(current_latitude, current_lon_diff)
        lon_diff = EarthEllipsoid.meters_to_parallel(current_latitude, parallel_length)
        recalculated_parallel_length = EarthEllipsoid.parallel_arc_length_with_sign(current_latitude, lon_diff)

        err_meters[k] = parallel_length - recalculated_parallel_length
        err_longitude[k] = current_lon_diff - lon_diff

    max_longitude_err = np.max(np.abs(np.array(err_longitude)))
    max_meters_err = np.max(np.abs(np.array(err_meters)))

    assert max_meters_err < meters_eps, "difference in meters is too big"
    assert max_longitude_err < geo_eps, "difference in radians is too big"


def test_geocentric_radius():
    # Checking geodetic_string_to_radians and geocentric_radius
    str_latitude = "0432635"
    str_longitude = "0395511"

    meters_eps = 1e-3

    latitude, longitude = EarthEllipsoid.geodetic_string_to_radians(str_latitude, str_longitude)
    geo_radius = EarthEllipsoid.geocentric_radius(latitude)  # in meters

    err = abs(6368070.2181 - geo_radius)

    assert err < meters_eps
