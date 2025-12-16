import numpy as np
import pytest

from CoordsStereo import StereoProjection
import CoordsCartesian3d


@pytest.fixture(scope="session")
def system_center():
    return [0.7582, 0.6967, 80.0]


@pytest.fixture(scope="session")
def get_point_geo_positions():

    number_of_points = 100

    # generate first and end points
    latitude1 = -np.pi / 2
    latitude2 = np.pi / 2

    longitude1 = -np.pi
    longitude2 = np.pi

    height1 = 0
    height2 = 13000

    # generate steps
    step_h = (height2 - height1) / number_of_points
    step_lat = (latitude2 - latitude1) / number_of_points
    step_lon = (longitude2 - longitude1) / number_of_points

    # generate points
    point_geo_pos = []
    for k in range(number_of_points):
        point_geo_pos.append([latitude1 + step_lat * k, longitude1 + step_lon * k, height1 + step_h * k])

    return point_geo_pos


@pytest.fixture(scope="session")
def get_point_geo_positions_around_the_center(system_center):
    center_geo_pos = system_center
    center_geometry = CoordsCartesian3d.PointGeometry(center_geo_pos)
    stereo_class = StereoProjection(center_geometry)

    number_of_points = 100
    max_range = 1000000

    # generate first and end points
    azimuth1 = 0
    azimuth2 = 2 * np.pi

    range1 = 100
    range2 = max_range

    height1 = center_geo_pos[2]
    height2 = 13000

    # generate steps
    step_h = (height2 - height1) / number_of_points
    step_az = (azimuth2 - azimuth1) / number_of_points
    step_range = (range2 - range1) / number_of_points

    # generate points
    points_geo_pos = []
    for k in range(number_of_points):
        r = range1 + step_range * k
        a = azimuth1 + step_az * k
        x = r * np.sin(a)
        y = r * np.cos(a)

        point_stereo_pos = [x, y, 0]
        point_geo_pos = stereo_class.stereo_to_geodetic(point_stereo_pos)
        point_geo_pos[2] = height1 + step_h * k
        points_geo_pos.append(point_geo_pos)

    return points_geo_pos