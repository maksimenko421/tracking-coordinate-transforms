import numpy as np
import sys
import time

import drawing_coords
import CoordsCartesian3d
from CoordsCartesianGeodesy import CartesianGeodesyCoordinates
from conftest import get_point_geo_positions_around_the_center


def test_geodetic_to_cartesian_local_and_back(system_center, get_point_geo_positions_around_the_center):
    # checking conversion from geographic coordinates to local cartesian coordinates and back.
    # Checking geodetic_to_cartesian_local and cartesian_local_to_geodetic functions

    graph_on = False  # switch on drawing

    geo_eps = 1e-10
    meters_eps = 1e-3

    center_geo_pos = system_center

    points_geo_pos = get_point_geo_positions_around_the_center
    number_of_points = len(points_geo_pos)


    err_geo = np.zeros([number_of_points, 2])
    err_local = np.zeros([number_of_points, 2])
    dist_from_center = np.zeros([number_of_points, ])

    center_geometry = CoordsCartesian3d.PointGeometry(center_geo_pos)
    cart_local_class = CartesianGeodesyCoordinates(center_geometry)

    run_time_geod_to_cart_local = 0
    run_time_cart_local_to_geod = 0

    for k in range(number_of_points):

        current_geo_pos = points_geo_pos[k][:]

        start_time_1 = time.perf_counter()
        point_cart_local_pos = cart_local_class.geodetic_to_local_cartesian(current_geo_pos)
        end_time_1 = time.perf_counter()
        run_time_geod_to_cart_local += end_time_1 - start_time_1

        dist_from_center[k] = np.sqrt(point_cart_local_pos[0] ** 2 + point_cart_local_pos[1] ** 2)

        start_time_2 = time.perf_counter()
        point_geo_pos1 = cart_local_class.cartesian_local_to_geodetic(point_cart_local_pos)
        end_time_2 = time.perf_counter()
        run_time_cart_local_to_geod += end_time_2 - start_time_2

        point_cart_local_pos1 = cart_local_class.geodetic_to_local_cartesian(point_geo_pos1)

        err_geo[k, :] = np.abs(np.array(current_geo_pos[:2]) - np.array(point_geo_pos1))
        err_local[k, :] = np.abs(np.array(point_cart_local_pos) - np.array(point_cart_local_pos1))

    max_lat_err = np.max(err_geo[:, 0])
    max_lon_err = np.max(err_geo[:, 1])

    max_x_err = np.max(err_local[:,0])
    max_y_err = np.max(err_local[:, 1])

    drawing_coords.plot_coord_cart_local_errors(graph_on, dist_from_center, err_geo, err_local)

    sys.stdout.write("\n")
    sys.stdout.write(">>> PRINT TEST RESULTS<<<\n")
    sys.stdout.write(f"run_time_geod_to_cart_local = {run_time_geod_to_cart_local}\n")
    sys.stdout.write(f"run_time_cart_local_to_geod = {run_time_cart_local_to_geod}\n")
    sys.stdout.flush()

    assert max_lat_err < geo_eps, "difference in latitude is too big"
    assert max_lon_err < geo_eps, "difference in longitude is too big"
    assert max_x_err < meters_eps, "difference in x-coordinate is too big"
    assert max_y_err < meters_eps, "difference in y-coordinate is too big"
