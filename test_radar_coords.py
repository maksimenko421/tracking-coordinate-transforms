import numpy as np
import sys
import time

import geometry_in_space
import drawing_coords
import CoordsCartesian3d
from CoordsRadar import RadarCoordinates


def test_geodetic_to_radar_and_back(system_center, get_point_geo_positions_around_the_center):
    # checking conversion from geographic coordinates to radar coordinates and back.
    # Checking geodetic_to_radar and radar_to_geodetic functions
    # max range = 1000km

    graph_on = False  # switch on drawing

    geo_eps = 1e-6
    azimuth_eps = 1e-5
    range_eps = 1e-6
    height_eps = 1e-7

    radar_geo_pos = system_center

    points_geo_pos = get_point_geo_positions_around_the_center
    number_of_points = len(points_geo_pos)

    err_geo = np.zeros([number_of_points, 3])
    err_polar = np.zeros([number_of_points, 3])
    center_geometry = CoordsCartesian3d.PointGeometry(radar_geo_pos)
    radar_class = RadarCoordinates(center_geometry)
    slant_range = np.zeros([number_of_points, ])

    run_time_radar_to_geod = 0
    run_time_geod_to_radar = 0
    for k in range(number_of_points):

        current_geo_pos = points_geo_pos[k][:]
        current_geo_pos[1] = geometry_in_space.normalization_to_main_period(current_geo_pos[1])

        point_polar_pos = radar_class.geodetic_to_radar(current_geo_pos)

        slant_range[k] = point_polar_pos[1]

        start_time_1 = time.perf_counter()
        point_geo_pos1 = radar_class.radar_to_geodetic(point_polar_pos)
        end_time_1 = time.perf_counter()
        run_time_radar_to_geod += end_time_1 - start_time_1

        start_time_2 = time.perf_counter()
        point_polar_pos1 = radar_class.geodetic_to_radar(point_geo_pos1)
        end_time_2 = time.perf_counter()
        run_time_geod_to_radar += end_time_2 - start_time_2

        err_geo[k, :] = np.abs(np.array(current_geo_pos) - np.array(point_geo_pos1))
        err_polar[k, :] = np.abs(np.array(point_polar_pos) - np.array(point_polar_pos1))
        if err_polar[k, 0] > np.pi:
            err_polar[k, 0] = np.abs(err_polar[k, 0] - 2 * np.pi)
        if err_polar[k, 0] < -np.pi:
            err_polar[k, 0] = np.abs(err_polar[k, 0] + 2 * np.pi)

    max_lat_err = np.max(err_geo[:, 0])
    max_lon_err = np.max(err_geo[:, 1])
    max_geo_height_err =  np.max(err_geo[:, 2])

    max_azimuth_err = np.max(err_polar[:,0])
    max_range_err = np.max(err_polar[:, 1])
    max_local_height_err = np.max(err_polar[:, 2])

    drawing_coords.plot_coord_radar_errors(graph_on, slant_range, err_geo, err_polar)

    sys.stdout.write(">>> PRINT TEST RESULTS<<<\n")
    sys.stdout.write(f"run_time_radar_to_geod = {run_time_radar_to_geod}\n")
    sys.stdout.write(f"run_time_geod_to_radar = {run_time_geod_to_radar}\n")
    sys.stdout.flush()

    assert max_lat_err < geo_eps, "difference in latitude is too big"
    assert max_lon_err < geo_eps, "difference in longitude is too big"
    assert max_azimuth_err < azimuth_eps, "difference in azimuth is too big"
    assert max_range_err < range_eps, "difference in slant range is too big"
    assert max_geo_height_err < height_eps, "difference in geo height is too big"
    assert max_local_height_err < height_eps, "difference in local height is too big"
