import numpy as np
import sys
import time

import drawing_coords
from CoordScaledGeodetic import ScaledGeodeticCoordinates


def test_geodetic_to_scaled_and_back(get_point_geo_positions, system_center):
    # checking conversion from geographic coordinates to local cartesian coordinates and back.
    # Checking geodetic_to_cartesian_local and cartesian_local_to_geodetic functions

    graph_on = False  # switch on drawing

    geo_eps = 1e-10
    meters_eps = 1e-8
    height_eps = 1e-12

    points_geo_pos = get_point_geo_positions
    number_of_points = len(points_geo_pos)

    err_geo = np.zeros([number_of_points, 3])
    err_local = np.zeros([number_of_points, 3])

    scaled_geodetic_class = ScaledGeodeticCoordinates(system_center)

    run_time_geod_to_sc_geod = 0
    run_time_sc_geod_to_geod = 0

    for k in range(number_of_points):

        current_geo_pos = points_geo_pos[k][:]

        start_time_1 = time.perf_counter()
        point_scaled_pos = scaled_geodetic_class.geodetic_to_scaled_geodetic(current_geo_pos)
        end_time_1 = time.perf_counter()
        run_time_geod_to_sc_geod += end_time_1 - start_time_1

        start_time_2 = time.perf_counter()
        point_geo_pos1 = scaled_geodetic_class.scaled_geodetic_to_geodetic(point_scaled_pos)
        end_time_2 = time.perf_counter()
        run_time_sc_geod_to_geod += end_time_2 - start_time_2

        point_scaled_pos1 = scaled_geodetic_class.geodetic_to_scaled_geodetic(point_geo_pos1)

        err_geo[k, :] = np.abs(np.array(current_geo_pos) - np.array(point_geo_pos1))
        err_local[k, :] = np.abs(np.array(point_scaled_pos) - np.array(point_scaled_pos1))

    max_lat_err = np.max(err_geo[:, 0])
    max_lon_err = np.max(err_geo[:, 1])
    max_geo_height_err =  np.max(err_geo[:, 2])

    max_x_err = np.max(err_local[:,0])
    max_y_err = np.max(err_local[:, 1])
    max_local_height_err = np.max(err_local[:, 2])

    drawing_coords.plot_coord_scaled_geo_errors(graph_on, points_geo_pos, err_geo, err_local)

    sys.stdout.write("\n")
    sys.stdout.write(">>> PRINT TEST RESULTS<<<\n")
    sys.stdout.write(f"run time_geod_to_sc_geod = {run_time_geod_to_sc_geod}\n")
    sys.stdout.write(f"run time_sc_geod_to_geod = {run_time_sc_geod_to_geod}\n")
    sys.stdout.flush()

    assert max_lat_err < geo_eps, "difference in latitude is too big"
    assert max_lon_err < geo_eps, "difference in longitude is too big"
    assert max_x_err < meters_eps, "difference in x-coordinate is too big"
    assert max_y_err < meters_eps, "difference in y-coordinate is too big"
    assert max_geo_height_err < height_eps, "difference in geo height is too big"
    assert max_local_height_err < height_eps, "difference in local height is too big"
