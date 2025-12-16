import numpy as np

from geodesy import EarthEllipsoid


def rotate_around_x(vector, phi):
    # rotation of a vector by an angle phi around the axis 0x
    vector = np.array(vector)
    sin_phi = np.sin(phi)
    cos_phi = np.cos(phi)

    rotation_matrix = np.array([[1, 0, 0],[0, cos_phi, -sin_phi],[0, sin_phi, cos_phi]])
    rotated_vector = rotation_matrix @ vector

    return rotated_vector


def rotate_around_y(vector, phi):
    # rotation of a vector by an angle phi around the axis 0y
    vector = np.array(vector)
    sin_phi = np.sin(phi)
    cos_phi = np.cos(phi)

    rotation_matrix = np.array([[cos_phi, 0, sin_phi], [0, 1, 0], [-sin_phi, 0, cos_phi]])
    rotated_vector = rotation_matrix @ vector

    return rotated_vector


def rotate_around_z(vector, phi):
    # rotation of a vector by an angle phi around the axis 0z
    vector = np.array(vector)
    sin_phi = np.sin(phi)
    cos_phi = np.cos(phi)

    rotation_matrix = np.array([[cos_phi, -sin_phi, 0],[sin_phi, cos_phi, 0], [0, 0, 1]])
    rotated_vector = rotation_matrix @ vector

    return rotated_vector


def get_rotation_matrix(vector, alpha):
    # finding the rotation matrix around the vector axis by the angle alpha counterclockwise

    co = np.cos(alpha)
    si = np.sin(alpha)

    vector = np.array(vector)
    norm_vector = np.linalg.norm(vector)
    x, y, z = vector / norm_vector

    matrix = [[co + (1 - co) * x * x, (1 - co) * x * y - si * z, (1 - co) * x * z + si * y],
              [(1 - co) * x * y + si * z, co + (1 - co) * y * y, (1 - co) * y * z - si * x],
              [(1 - co) * x * z - si * y, (1 - co) * y * z + si * x, co + (1 - co) * z * z]]

    return np.array(matrix)


def get_rotation_matrix2(axis, alpha):
    co = np.cos(alpha)
    si = np.sin(alpha)

    axis /= np.linalg.norm(axis)
    x, y, z = axis  # Assuming axis is already a unit vector

    K = np.array([[0, -z, y], [z, 0, -x], [-y, x, 0]])  # Skew-symmetric matrix

    matrix = np.eye(3) + si * K + (1 - co) * (K @ K)  # Rodrigues' formula

    return np.array(matrix)


def normal_to_ellipsoid(point_xyz_proj):
    # Input: point_xyz_proj is a Cartesian 3D coordinates of a point which belongs to the ellipsoid EarthEllipsoid
    # ellipsoid has equation x^2/a^2 + y^2/a^2 + z^2/b^2 = 1
    # direction of normal vector at given point point_xyz_proj coincides
    # with direction of partial derivative vector to ellipsoid at given point

    a2 = EarthEllipsoid.a2  # major semi-axis squared
    b2 = EarthEllipsoid.b2  # minor semi-axis squared

    x0, y0, z0 = point_xyz_proj  # the point must belong to the ellipsoid

    n0 = [x0, y0, a2 * z0 / b2]  # normal to ellipsoid at point
    n0 = n0 / np.linalg.norm(n0)  # normalization

    return n0


def vector_towards_the_nearest_pole(point_xyz_proj,  zero_eps = 1e-14):
    # Input: point_xyz_pos = [x0,y0,z0] is Cartesian 3D coordinates of the point at the ellipsoid
    # Output: 3D unit vector of direction towards the North or South Pole in the tangent plane to the EarthEllipsoid
    # the ellipsoid has the equation x^2/a^2 + y^2/a^2 + z^2/b^2 = 1
    # equation of the tangent plane at the point is
    # x * x0 + y * y0 + z * z0 * a^2 / b^2 = a^2
    # finding the intersection point of the tangent plane and the 0z axis (substituting x=0, y=0),
    # we get z = b^2 / z0, so, the intersection point is (0, 0, b^2 / z0)
    # to get the coordinates of the vector, subtract from (0, 0, b^2 / z0) the vector (x0, y0, z0)

    b2 = EarthEllipsoid.b2  # square of the earth's minor semi-axis in meters
    x0, y0, z0 = point_xyz_proj

    if np.abs(z0) > zero_eps:  # the point is far enough from the Equator
        q0 = [-x0, -y0, (b2 / z0 - z0)]
        q0 = q0 / np.linalg.norm(q0)
    else:
        q0 = np.array([0,0,1]) if z0 >= 0 else np.array([0,0,-1])

    return q0


def normal_to_the_meridian_plane(point_xyz_proj):
    # Input: point_xyz_proj should belong to ellipsoid
    # Output: unit vector of the normal to the meridian section.
    # Meridian section is the plane containing the Earth center, Poles and point_xyz_proj

    q0 = vector_towards_the_nearest_pole(point_xyz_proj)
    n = normal_to_ellipsoid(point_xyz_proj)

    unit_normal = np.cross(q0, n) # normal to the desired plane

    return unit_normal


def normal_plane_at_azimuth_direction(point_xyz_proj, azimuth):
    # Input: point_xyz_proj must belong to the ellipsoid
    # azimuth is the angle in which the direction the desired plane will be
    # Output: coefficients of the desired plane

    x0, y0, z0 = point_xyz_proj

    # Step 1. Finding meridian plane at the radar projection point
    meridian_normal = normal_to_the_meridian_plane(point_xyz_proj)

    # Step 2. Finding rotation matrix to turn around normal to ellipsoid by point_azimuth angle
    unit_normal = normal_to_ellipsoid(point_xyz_proj)
    rotation_matrix = get_rotation_matrix(unit_normal, -azimuth)

    # Step 3. Finding the plane of normal section of  the ellipsoid at the radar projection point which contains points M and A (A is the projection of M)
    sx, sy, sz = rotation_matrix @ meridian_normal
    s0 = -(sx * x0 + sy * y0 + sz * z0)  # free coefficient of the plane

    return [sx, sy, sz, s0]


def check_point_on_ellipsoid(point_xyz_pos):
    a2 = EarthEllipsoid.a2
    b2 = EarthEllipsoid.b2
    x, y, z = point_xyz_pos
    res = x * x / a2 + y * y / a2 + z * z / b2 - 1
    return res


def check_point_on_plane(point_xyz_pos, plane_coefficients):
    x, y, z = point_xyz_pos
    sx, sy, sz, s0 = plane_coefficients
    res = sx * x + sy * y + sz * z + s0
    return res


def project_onto_ellipsoid_Lagrange(point_xyz_pos, zero_eps=1e-10, max_iter=20):
    # project onto the ellipsoid using Lagrange Multipliers method
    a2 = EarthEllipsoid.a2
    b2 = EarthEllipsoid.b2
    x, y, z = point_xyz_pos  # initial guess

    for _ in range(max_iter):
        fx = x**2 / a2 + y**2 / a2 + z**2 / b2 - 1
        dfx = 2 * (x / a2 + y / a2 + z / b2)
        lambda_ = fx / dfx
        x -= lambda_ * x / a2
        y -= lambda_ * y / a2
        z -= lambda_ * z / b2
        if abs(fx) < zero_eps:
            break
    return np.array([x, y, z])


def project_onto_ellipsoid(point_xyz_pos, zero_eps=1e-1):
    # Finding the point (xe, ye, ze) - the projection of the point point_xyz_pos = (xa, ya, za) onto the ellipsoid
    # as the intersection point of the line passing through the point and perpendicular to the ellipsoid,
    # and the surface of the ellipsoid, by solving the 4th degree equation
    # with respect to xe with coefficients poly_coefs = [1, tmp2, tmp3, tmp4, tmp5];
    # the most precisely method

    a2 = EarthEllipsoid.a2  # major semi-axis squared in meters
    b = EarthEllipsoid.b  # minor semi-axis in meters
    e2 = EarthEllipsoid.e2  # eccentricity squared
    e4 = e2 * e2

    xa, ya, za = point_xyz_pos
    r2 = xa * xa + ya * ya  # auxiliary constant

    ## case 1: |xa| > 0
    if np.abs(xa) > zero_eps:
        coef3 = -2 * xa / e2  # coefficient of the 3rd power
        coef2 = xa * xa * (r2 + (1 - e2) * za * za - a2 * e4) / (e4 * r2)  # coefficient of the 2nd power
        coef1 = 2 * a2 * xa ** 3 / (e2 * r2)  # coefficient of the 1st power
        coef0 = -a2 * xa ** 4 / (e4 * r2)  # free coefficient

        # finding polynomial roots
        poly_coefs = [1, coef3, coef2, coef1, coef0]
        poly_roots = np.roots(poly_coefs)

        # discard unnecessary solutions: complex solutions are not suitable and
        # the intersection point must have the same sign as number xa
        # moreover we need the line which intersects polynomial in 2 points, therefore
        # the polynomial root with opposite sign must exist
        xe = 0
        for poly_root in poly_roots:
            if np.abs(poly_root.imag) < zero_eps and xa * poly_root.real > 0 and poly_root.real < 1e+7:
                # if any(np.abs(poly_root.real + s.real) < 10 for s in poly_roots):
                xe = poly_root.real  # suitable solution
                break

        ye = ya * xe / xa
        ze = za * (1 - e2) * xe / (xa - e2 * xe)

        proj_point = [xe, ye, ze]

    elif np.abs(ya) > zero_eps:

        coef3 = -2 * ya / e2  # coefficient of the 3rd power
        coef2 = ya * ya * (r2 + (1 - e2) * za * za - a2 * e4) / (e4 * r2)  # coefficient of the 2nd power
        coef1 = 2 * a2 * ya ** 3 / (e2 * r2)  # coefficient of the 1st power
        coef0 = -a2 * ya ** 4 / (e4 * r2)  # free coefficient

        # finding polynomial roots
        poly_coefs = [1, coef3, coef2, coef1, coef0]
        poly_roots = np.roots(poly_coefs)

        # discard unnecessary solutions: complex solutions are not suitable and
        # the intersection point must have the same sign as number xa
        # moreover we need the line which intersects polynomial in 2 points, therefore
        # the polynomial root with opposite sign must exist
        ye = 0
        for poly_root in poly_roots:
            if np.abs(poly_root.imag) < zero_eps and ya * poly_root.real > 0 and poly_root.real < 1e+7:
                # if any(np.abs(poly_root.real + s.real) < 10 for s in poly_roots):
                ye = poly_root.real  # suitable solution
                break

        xe = xa * ye / ya
        ze = za * (1 - e2) * ye / (ya - e2 * ye)

        proj_point = [xe, ye, ze]

        """
        coef3 = -2 * ya / e2  # coefficient of the 3rd power
        coef2 = ya * ya * (ya * ya + (1 - e2) * za * za - a2 * e4) / e4  # coefficient of the 2nd power
        coef1 = 2 * a2 * ya / e2  # coefficient of the 1st power
        coef0 = -a2 * ya ** 2 / e4  # free coefficient

        # finding polynomial roots
        poly_coefs = [1, coef3, coef2, coef1, coef0]
        poly_roots = np.roots(poly_coefs)

        # discard unnecessary solutions: complex solutions are not suitable and
        # the intersection point must have the same sign as number xa
        # moreover we need the line which intersects polynomial in 2 points, therefore
        # the polynomial root with opposite sign must exist
        ye = 0
        for poly_root in poly_roots:
            if np.abs(poly_root.imag) < zero_eps and ya * poly_root.real > 0 and poly_root.real < 1e+7:
                # if any(np.abs(poly_root.real + s.real) < 10 for s in poly_roots):
                ye = poly_root.real  # suitable solution
                break

        ze = za * (1 - e2) * ye / (ya - e2 * ye)
        proj_point = [0, ye, ze]
    """
    else:  # xa = ya = 0, one of the poles
        ze = np.sign(za) * b
        proj_point = [0, 0, ze]

    return proj_point


def get_rotation_matrix_for_enu(center_geo_pos):

    lat = center_geo_pos[0]
    lon = center_geo_pos[1]

    cos_lat = np.cos(lat)
    sin_lat = np.sin(lat)
    cos_lon = np.cos(lon)
    sin_lon = np.sin(lon)

    matrix = [[-sin_lon, cos_lon, 0],
              [-sin_lat * cos_lon, -sin_lat * sin_lon, cos_lat],
              [cos_lat * cos_lon, cos_lat * sin_lon, sin_lat]]

    return np.array(matrix)


def normalization_to_main_period(x):
    # normalization to the interval [0, 2 * pi)

    x = x % (2 * np.pi)  # shifting x to [0, 2*pi) for positive or to (-2*pi,0] for negative numbers

    # shifting result to [0, 2 pi)
    if x < 0:
        x += 2 * np.pi

    return x
