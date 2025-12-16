import numpy as np

import geometry_in_space
import drawing_coords
from CoordsCartesian3d import Cartesian3d


def test_geodetic_to_cartesian3d_and_back(get_point_geo_positions):
    # checking conversion from Geodetic coordinates to Cartesian 3D coordinates and back.
    # Checking geodetic_to_cartesian3D and cartesian3D_to_geodetic methods

    graph_on = False  # switch on drawing

    geo_eps = 1e-14
    meters_eps = 1e-7
    height_eps = 1e-7

    point_geo_pos = get_point_geo_positions
    number_of_points = len(point_geo_pos)

    err_geo = np.zeros([number_of_points, 3])
    err_cart3d = np.zeros([number_of_points, 3])

    cartesian3d_class = Cartesian3d()

    for k in range(number_of_points):

        current_geo_pos = point_geo_pos[k][:]
        point_xyz_pos = cartesian3d_class.geodetic_to_cartesian3d(current_geo_pos)
        point_geo_pos1 = cartesian3d_class.cartesian3d_to_geodetic(point_xyz_pos)
        point_xyz_pos1 = cartesian3d_class.geodetic_to_cartesian3d(point_geo_pos1)

        current_geo_pos[1] = geometry_in_space.normalization_to_main_period(current_geo_pos[1])
        err_geo[k, :] = np.abs(np.array(current_geo_pos) - np.array(point_geo_pos1))
        err_cart3d[k, :] = np.abs(np.array(point_xyz_pos) - np.array(point_xyz_pos1))

    max_lat_err = np.max(err_geo[:, 0])
    max_lon_err = np.max(err_geo[:, 1])
    max_geo_height_err = np.max(err_geo[:, 2])

    max_x_err = np.max(err_cart3d[:, 0])
    max_y_err = np.max(err_cart3d[:, 1])
    max_z_err = np.max(err_cart3d[:, 2])


    drawing_coords.plot_coord_cart3d_errors(graph_on, point_geo_pos, err_geo, err_cart3d)


    assert max_x_err < meters_eps, "difference in x-coordinate is too big"
    assert max_y_err < meters_eps, "difference in y-coordinate is too big"
    assert max_z_err < meters_eps, "difference in z-coordinate is too big"
    assert max_lat_err < geo_eps, "difference in latitude is too big"
    assert max_lon_err < geo_eps, "difference in longitude is too big"
    assert max_geo_height_err < height_eps, "difference in geo height is too big"
