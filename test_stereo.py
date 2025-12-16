import numpy as np
import sys
import time

import geometry_in_space
import drawing_coords
import CoordsCartesian3d
from CoordsStereo import StereoProjection


def test_geodetic_to_stereo_and_back(system_center, get_point_geo_positions_around_the_center):
    # checking conversion from geographic coordinates to stereo coordinates and back.
    # Checking geodetic_to_stereo and stereo_to_geodetic functions
    # max range = 1000km

    graph_on = False  # switch on drawing

    geo_eps = 1e-9
    meters_eps = 1e-5
    height_eps = 1e-12

    radar_geo_pos = system_center

    points_geo_pos = get_point_geo_positions_around_the_center
    number_of_points = len(points_geo_pos)

    err_geo = np.zeros([number_of_points, 3])
    err_stereo = np.zeros([number_of_points, 3])
    stereo_range = np.zeros([number_of_points, ])

    center_geometry = CoordsCartesian3d.PointGeometry(radar_geo_pos)
    stereo_class = StereoProjection(center_geometry)

    run_time_geod_to_stereo = 0
    run_time_stereo_to_geod = 0
    for k in range(number_of_points):

        current_geo_pos = points_geo_pos[k][:]

        start_time_1 = time.perf_counter()
        point_stereo_pos = stereo_class.geodetic_to_stereo(current_geo_pos)
        end_time_1 = time.perf_counter()
        run_time_geod_to_stereo += end_time_1 - start_time_1

        stereo_range[k] = np.sqrt(point_stereo_pos[0] ** 2 + point_stereo_pos[1] ** 2 )

        start_time_2 = time.perf_counter()
        point_geo_pos1 = stereo_class.stereo_to_geodetic(point_stereo_pos)
        end_time_2 = time.perf_counter()
        run_time_stereo_to_geod += end_time_2 - start_time_2

        point_stereo_pos1 = stereo_class.geodetic_to_stereo(point_geo_pos1)

        current_geo_pos[1] = geometry_in_space.normalization_to_main_period(current_geo_pos[1])

        err_geo[k, :] = np.abs(np.array(current_geo_pos) - np.array(point_geo_pos1))
        err_stereo[k, :] = np.abs(np.array(point_stereo_pos) - np.array(point_stereo_pos1))

    max_lat_err = np.max(err_geo[:, 0])
    max_lon_err = np.max(err_geo[:, 1])
    max_geo_height_err =  np.max(err_geo[:, 2])

    max_x_err = np.max(err_stereo[:,0])
    max_y_err = np.max(err_stereo[:, 1])
    max_local_height_err = np.max(err_stereo[:, 2])

    drawing_coords.plot_coord_stereo_errors(graph_on, stereo_range, err_geo, err_stereo)

    sys.stdout.write("\n")
    sys.stdout.write(">>> PRINT TEST RESULTS<<<\n")
    sys.stdout.write(f"run_time_geod_to_stereo = {run_time_geod_to_stereo}\n")
    sys.stdout.write(f"run_time_stereo_to_geod = {run_time_stereo_to_geod}\n")
    sys.stdout.flush()

    assert max_lat_err < geo_eps, "difference in latitude is too big"
    assert max_lon_err < geo_eps, "difference in longitude is too big"
    assert max_x_err < meters_eps, "difference in x-coordinate is too big"
    assert max_y_err < meters_eps, "difference in y-coordinate is too big"
    assert max_geo_height_err < height_eps, "difference in geo height is too big"
    assert max_local_height_err < height_eps, "difference in local height is too big"
