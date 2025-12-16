import numpy as np

import geometry_in_space
from CoordsCartesian3d import Cartesian3d, PointGeometry
from CoordsRadar import RadarCoordinates


def test_normal_plane_containing_the_meridian():

    meters_eps = 1e-8

    point_geo_pos = [-1.55814, 0, 80.0]
    point_geo_proj = point_geo_pos
    point_geo_proj[2] = 0
    point_xyz_proj = Cartesian3d.geodetic_to_cartesian3d(point_geo_proj)

    normal_plane_1 = geometry_in_space.normal_to_the_meridian_plane(point_xyz_proj)

    test_normal_1 = [0,1,0]

    err_xyz_1 = np.abs(np.cross(np.array(normal_plane_1), np.array(test_normal_1)))

    point_geo_pos = [-1.55814, np.pi/2, 80.0]
    point_geo_proj = [point_geo_pos[0], point_geo_pos[1], 0]
    point_xyz_proj = Cartesian3d.geodetic_to_cartesian3d(point_geo_proj)

    normal_plane_2 = geometry_in_space.normal_to_the_meridian_plane(point_xyz_proj)

    test_normal_2 = [-1, 0, 0]

    err_xyz_2 = np.abs(np.cross(np.array(normal_plane_2) , np.array(test_normal_2)))

    max_err_1 = np.max(err_xyz_1)
    max_err_2 = np.max(err_xyz_2)

    assert max_err_1 < meters_eps, "difference in coordinates is too big in plane 1"
    assert max_err_2 < meters_eps, "difference in coordinates is too big in plane 2"


def test_vector_to_the_pole_and_normal():

    meters_eps = 1e-8

    point_geo_pos = [1.55814, 0.6967, 80.0]
    point_xyz_pos = Cartesian3d.geodetic_to_cartesian3d(point_geo_pos)
    point_geo_proj = [point_geo_pos[0], point_geo_pos[1], 0]

    point_xyz_proj = Cartesian3d.geodetic_to_cartesian3d(point_geo_proj)

    vector_to_Pole = geometry_in_space.vector_towards_the_nearest_pole(point_xyz_proj)
    normal = geometry_in_space.normal_to_ellipsoid(point_xyz_pos)

    inner_prod = np.dot(vector_to_Pole, normal)
    abs_err = np.abs(inner_prod)

    assert abs_err < meters_eps, "difference in coordinates is too big in plane 1"


def test_vector_rotation():

    meters_eps = 1e-8

    vector_a = np.array([1,2,3], dtype=float)
    vector_b = np.array([2,4,1], dtype=float)

    axis = np.cross(vector_a, vector_b)
    axis_norm = np.linalg.norm(axis)

    # Handle collinear vectors (rotation is undefined in this case)
    if axis_norm < 1e-12:
        raise ValueError("Rotation axis is ill-defined: input vectors might be collinear.")
    axis /= axis_norm

    cos_phi = np.dot(vector_a, vector_b) / (np.linalg.norm(vector_a) * np.linalg.norm(vector_b))
    cos_phi = np.clip(cos_phi, -1.0, 1.0)  # Avoid floating-point issues
    phi = np.arccos(cos_phi)

    # Compute the rotation matrix and apply it
    rotation_matrix = np.array(geometry_in_space.get_rotation_matrix(axis, phi))
    det_rm = np.linalg.det(rotation_matrix)
    rotation_matrix2 = np.array(geometry_in_space.get_rotation_matrix2(axis, phi))
    vector_b1 = np.matmul(rotation_matrix, vector_a) * np.linalg.norm(vector_b) / np.linalg.norm(vector_a)

    # Compute the error
    err_xyz_b = np.linalg.norm(vector_b1 - vector_b)

    assert err_xyz_b < meters_eps, "difference in coordinates is too big"


def test_specific_rotations():
    meters_eps = 1e-8

    vector_a = np.array([1, 1, 1], dtype=float)
    phi = np.pi / 4

    vector_x = geometry_in_space.rotate_around_x(vector_a, phi)
    vector_y = geometry_in_space.rotate_around_y(vector_a, phi)
    vector_z = geometry_in_space.rotate_around_z(vector_a, phi)

    err_x = np.linalg.norm(vector_x - np.array([1, 0, np.sqrt(2)]))
    err_y = np.linalg.norm(vector_y - np.array([np.sqrt(2), 1, 0]))
    err_z = np.linalg.norm(vector_z - np.array([0, np.sqrt(2), 1]))

    assert err_x < meters_eps, "difference in coordinates is too big"
    assert err_y < meters_eps, "difference in coordinates is too big"
    assert err_z < meters_eps, "difference in coordinates is too big"


def test_normal_plane_at_azimuth_direction():
    meters_eps = 1e-10

    radar_geo_pos = [1.55814, 0.6967, 0.0]
    point_geo_pos = [1.48, 0.66892, 12192.0]  # distance 500 km

    #radar_geo_pos = [1.55814, 0.6967, 200.0]
    radar_geo_proj = [radar_geo_pos[0], radar_geo_pos[1], 0]
    radar_xyz_proj = Cartesian3d.geodetic_to_cartesian3d(radar_geo_proj)

    #point_geo_pos = [1.48, 0.66892, 12192.0]  # distance 500 km
    x0, y0, z0 = Cartesian3d.geodetic_to_cartesian3d(point_geo_pos)

    radar_geometry = PointGeometry(radar_geo_pos)
    radar_class = RadarCoordinates(radar_geometry)
    point_polar_pos = radar_class.geodetic_to_radar(point_geo_pos)
    point_azimuth = point_polar_pos[0]

    # function call
    [sx, sy, sz, s0] = geometry_in_space.normal_plane_at_azimuth_direction(radar_xyz_proj, point_azimuth)

    # checking that point x0, y0, z0 belongs to the plane
    err_xyz = np.abs(sx * x0 + sy * y0 + sz * z0 + s0)

    assert err_xyz < meters_eps, "difference in coordinates is too big"


def test_orthogonality():
    meters_eps = 1e-10

    radar_geo_pos = [1.55814, 0.6967, 1200.0]
    radar_geo_proj = [radar_geo_pos[0], radar_geo_pos[1], 0]
    radar_xyz_proj = Cartesian3d.geodetic_to_cartesian3d(radar_geo_proj)

    # function call
    n = geometry_in_space.normal_to_ellipsoid(radar_xyz_proj)
    q0 = geometry_in_space.vector_towards_the_nearest_pole(radar_xyz_proj)

    # checking
    err_xyz = np.abs(np.dot(n, q0))

    assert err_xyz < meters_eps, "difference in coordinates is too big"
