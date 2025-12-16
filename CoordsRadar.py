import numpy as np

from geodesy import EarthEllipsoid
from CoordsCartesian3d import Cartesian3d
import geometry_in_space

# Radar coordinates: [azimuth, range, altitude] relatively radar position in [radians, meters, meters]
# azimuth is the angle in the tangent plane between the direction to a point and the direction to the north, measured clockwise
# range is the distance between the radar and the point
# altitude is the distance from the point to the surface of the ellipsoid

# Geodetic coordinates [latitude, longitude, altitude] in [radians, radians, meters]
# the latitude is an acute angle formed by the normal to the surface of an ellipsoid and the plane of the equator,
# belongs to [0, pi/2] in the northern hemisphere and to [-pi/2,0] in the southern hemisphere.
# the longitude is a dihedral angle formed by the plane of the prime meridian and a given meridian,
# belongs to [0, pi] in the eastern hemisphere and to [-pi,0] in the western hemisphere.
# the altitude is the distance from the point to the surface of the ellipsoid


class RadarCoordinates:
    def __init__(self, radar_geometry):
        self.point_geo_pos = np.zeros([3, ])
        self.point_polar_pos = np.zeros([3, ])
        self.radar_geo_pos = radar_geometry.geo_pos
        self.radar_geometry = radar_geometry


    def radar_to_geodetic_sphere(self, point_polar_pos, radius, zero_eps=1e-4):
        # Input: radar_geo_pos [latitude, longitude, altitude] are Geodetic coordinates of the radar position
        # point_polar_pos[azimuth, range, altitude] are radar coordinates of the point
        # radius - geocentric radius of the sphere in the radar position
        # Output: point_geo_pos[latitude, longitude, altitude] are Geodetic coordinates of the point on a sphere in the
        # northern hemisphere

        radar_geo_pos = self.radar_geo_pos

        # calculating in advance what will be useful
        sin_radar_latitude = np.sin(radar_geo_pos[0])
        cos_radar_latitude = np.cos(radar_geo_pos[0])
        radar_longitude = radar_geo_pos[1]
        azimuth, distance, altitude = point_polar_pos

        # Step 1. Applying the law of cosine to a triangle with vertices: center of the Earth, radar, point
        c = altitude + radius  # distance from the center of the Earth to the point
        d = radius + radar_geo_pos[2]  # distance from the center of the earth to the radar
        cos_phi0 = (d * d + c * c - distance * distance) / (2 * d * c)  # phi0 - distance between radar and point on unit sphere
        if np.abs(cos_phi0) > 1:
            print('radar_to_geodetic_sphere: invalid point')
            return None, None
        sin_phi0 = np.sqrt(1 - cos_phi0 * cos_phi0)

        # Step 2. Applying spherical law of cosines for a spherical triangle on the unit sphere with vertices: North Pole, radar, point
        sin_latitude = sin_radar_latitude * cos_phi0 + cos_radar_latitude * sin_phi0 * np.cos(azimuth)
        cos_latitude = np.sqrt(1 - sin_latitude * sin_latitude)
        latitude = np.asin(sin_latitude)

        if not np.isnan(latitude):
            if abs(cos_latitude) > zero_eps:  # the point far enough from the pole
                longitude = radar_longitude + np.asin(sin_phi0 * np.sin(azimuth) / cos_latitude)  # spherical law of cosines
                longitude = geometry_in_space.normalization_to_main_period(longitude)
            else:  # handling pole crossing
                print('the point near the Pole')
                longitude = 0
        else:
            longitude = np.nan

        return [latitude, longitude, altitude]


    def radar_to_geodetic(self, point_polar_pos, max_iter=20, zero_eps=1e-10):
        # Input: radar_geometry contains also Geodetic coordinates of radar position
        # point_polar_pos = [azimuth, range, altitude] are polar coordinates of the point in the radar coordinate system
        # Output: point_geo_pos = [latitude, longitude, altitude] are Geodetic coordinates of the point

        # Step 0. Auxiliary
        a2 = EarthEllipsoid.a2
        b2 = EarthEllipsoid.b2
        a2_b2 = EarthEllipsoid.a2_b2
        a4_b4 = EarthEllipsoid.a4_b4
        radar_geometry = self.radar_geometry
        radius = radar_geometry.geocentric_radius
        xr, yr, zr = radar_geometry.xyz_pos

        point_azimuth, point_range, point_height = point_polar_pos
        point_range_sq = point_range * point_range
        point_height_sq = point_height * point_height

        # Step 1. Finding the normal plane at the azimuth direction at the radar projection point
        [sx, sy, sz, s0] = geometry_in_space.normal_plane_at_azimuth_direction(radar_geometry.xyz_proj, point_azimuth)

        # Step 2. A(xe, ye, ze) is the unknown projection point of the desired point M(x0, y0, z0) onto ellipsoid
        # Finding initial guess for the system relatively of coordinates of points M(x0,y0,z0) and A(xe,ye,ze).
        cartesian3d_class = Cartesian3d()
        point_geo_pos_sphere = self.radar_to_geodetic_sphere(point_polar_pos, radius)
        x00, y00, z00 = cartesian3d_class.geodetic_to_cartesian3d(point_geo_pos_sphere)
        point_polar_proj = [point_polar_pos[0], point_polar_pos[1], 0]
        point_geo_proj_sphere = self.radar_to_geodetic_sphere(point_polar_proj, radius)
        xe0, ye0, ze0 = cartesian3d_class.geodetic_to_cartesian3d(point_geo_proj_sphere)

        if np.abs(x00) > zero_eps:
            k0 = x00 / xe0 - 1  # x0 = (k+1) * xe
        elif np.abs(y00) > zero_eps:
            k0 = y00 / ye0 - 1  # y0 = (k+1) * ye
        else:
            k0 = (z00 / ze0 - 1) / a2_b2  # z0 = (k * a2_b2 + 1) * ze

        """
        initial_guess = np.array([xe0, ye0, ze0, k0])
        def system_of_equations(vars):
            xe, ye, ze, k = vars
            eq1 = sx * xe + sy * ye + sz * ze + s0  #  equation for point M azimuth
            eq2 = xe ** 2 / a2 + ye ** 2 / a2 + ze ** 2 / b2 - 1 # equation for point A belongs to the Earth ellipsoid
            eq3 = k ** 2 * (xe ** 2 + ye ** 2 + ze ** 2 * a4_b4) - point_height_sq  # equation for point M height
            eq4 = ((k + 1) * xe - xr) ** 2 + ((k + 1) * ye - yr) ** 2 + ((k * a2_b2 + 1) * ze - zr) ** 2 - point_range_sq  # equation for point M slant range
            return [eq1, eq2, eq3, eq4]
    
        x = optimize.fsolve(system_of_equations, initial_guess, epsfcn=zero_eps)
        """

        def F(x):
            # Define the system of quadratic equations.
            system_function = [
                sx * x[0] + sy * x[1] + sz * x[2] + s0,  # First equation
                x[0] ** 2 / a2 + x[1] ** 2 / a2 + x[2] ** 2 / b2 - 1,  # Second equation
                x[3] ** 2 * (x[0] ** 2 + x[1] ** 2 + x[2] ** 2 * a4_b4) - point_height_sq,  # equation for point M height
                ((x[3] + 1) * x[0] - xr) ** 2 + ((x[3] + 1) * x[1] - yr) ** 2 + (
                            (x[3] * a2_b2 + 1) * x[2] - zr) ** 2 - point_range_sq]  # equation for point M slant range
            return np.array(system_function)

        def J_F(x):
            # Compute the Jacobian matrix.

            j33 = 2 * ((x[3] + 1) * x[0] - xr) * x[0] + 2 * ((x[3] + 1) * x[1] - yr) * x[1] + 2 * (
                        (x[3] * a2_b2 + 1) * x[2] - zr) * x[2] * a2_b2
            system_jacobian = [[sx, sy, sz, 0],
                               [2 * x[0] / a2, 2 * x[1] / a2, 2 * x[2] / b2, 0],
                               [2 * x[0] * x[3] ** 2, 2 * x[1] * x[3] ** 2, 2 * x[2] * x[3] ** 2 * a4_b4,
                                2 * x[3] * (x[0] ** 2 + x[1] ** 2 + x[2] ** 2 * a4_b4)],
                               [2 * ((x[3] + 1) * x[0] - xr) * (x[3] + 1), 2 * ((x[3] + 1) * x[1] - yr) * (x[3] + 1),
                                2 * ((x[3] * a2_b2 + 1) * x[2] - zr) * (x[3] * a2_b2 + 1), j33]]

            return np.array(system_jacobian)

        # Newton Raphson method for solving quadratic system
        x = [xe0, ye0, ze0, k0]
        for iter in range(max_iter):
            J = J_F(x)  # Jacobian matrix
            F_val = F(x)  # Function values
            # delta_x0 = np.linalg.solve(J, -F_val)  # Solve the linear system J * delta_x = -F

            det_J = np.linalg.det(J)
            if np.abs(det_J) > zero_eps:
                inv_jacobian = np.linalg.inv(J)
                delta_x = inv_jacobian @ F_val
            else:
                delta_x = np.zeros(4)

            x -= delta_x  # Update estimate

            if np.linalg.norm(delta_x) < zero_eps:
                break

        xe, ye, ze, k = x
        x0 = (k + 1) * xe
        y0 = (k + 1) * ye
        z0 = (k * a2_b2 + 1) * ze

        point_geo_pos = cartesian3d_class.cartesian3d_to_geodetic([x0, y0, z0])
        point_geo_pos[2] = point_polar_pos[2]

        return point_geo_pos


    def geodetic_to_radar(self, point_geo_pos):
        # Input: radar_geo_pos = (latitude, longitude, altitude) - vector of radar Geodetic coordinates;
        # point_geo_pos = (latitude, longitude, altitude) - vector of Geodetic coordinates on the ellipsoid;
        # Output: point_polar_pos = (azimuth, range, altitude) vector coordinates relatively radar position
        # (azimuth is in radians, range, altitude are in meters)
        # The result is accurate. the function runs 4 times faster than geodetic_to_radar_auxiliary

        radar_geometry = self.radar_geometry
        radar_geo_pos = self.radar_geo_pos
        radar_latitude = radar_geo_pos[0]
        radar_longitude = radar_geo_pos[1]
        radar_xyz_pos = radar_geometry.xyz_pos  # radar Cartesian 3D coordinates
        radar_xyz_proj = radar_geometry.xyz_proj  # radar projection Cartesian 3D coordinates

        cartesian3d_class = Cartesian3d()
        point_xyz_pos = cartesian3d_class.geodetic_to_cartesian3d(point_geo_pos)  # point Cartesian 3D coordinates

        vector_to_point = np.array(point_xyz_pos) - np.array(
            radar_xyz_pos)  # vector of direction to the point from the radar point
        point_range = np.linalg.norm(vector_to_point)  # distance to point

        vector_from_radar_proj = np.array(point_xyz_pos) - np.array(
            radar_xyz_proj)  # vector of direction to the point from the radar projection point

        # alignment of the tangent plane at the radar projection position with the xy plane, that is parallel to the Equator
        # after this the normal vector will be directed along the z axis, and the vector to north in tangent plane will be directed along the x-axis,
        # and we calculate azimuth in 2D plane xOy
        point_vector = geometry_in_space.rotate_around_z(vector_from_radar_proj,
                                                         -radar_longitude)  # turning radar position to Greenwich longitude
        point_vector = geometry_in_space.rotate_around_y(point_vector,
                                                         radar_latitude - 0.5 * np.pi)  # turning tangent plane to be parallel to the Equator

        # after the space rotation vector to nord in tangent plane will be co-directed with the vector [-1, 0, 0]]
        point_azimuth = np.pi - np.atan2(point_vector[1], point_vector[0])  # belongs to [0, 2 * pi)
        point_azimuth = geometry_in_space.normalization_to_main_period(point_azimuth)

        point_polar_pos = [point_azimuth, point_range, point_geo_pos[2]]

        return point_polar_pos


    def radar_to_geodetic_covariance(self, point_polar_pos, cov_radar):
        # recalculation covariance matrix of the point from radar coordinates to geodetic coordinates
        # input:
        # radar_geometry is a structure which contains geometric parameters of the radar
        # point_polar_pos is a point position in radar coordinates:
        # azimuth (in radians), range (in meters), and altitude (in meters),
        # cov_radar is a 2x2 diagonal matrix of radar variances in azimuth and range,
        # output 2x2 covariance matrix in geodetic coordinates

        radar_geometry = self.radar_geometry
        radar_geo_pos = self.radar_geo_pos
        radar_cos_latitude = radar_geometry.cos_latitude
        radar_sin_latitude = radar_geometry.sin_latitude

        radius = radar_geometry.radius  # radius of the sphere at the radar point
        radius_squared = radius * radius

        azimuth = point_polar_pos[0]  # point azimuth
        sin_azimuth = np.sin(azimuth)
        cos_azimuth = np.cos(azimuth)
        range = point_polar_pos[1]  # point range
        range_squared = range * range

        # calculating derivatives to find the Jacobian matrix ja
        point_radius = radius + point_polar_pos[2] - radar_geo_pos[2]  # radius to point
        point_radius_squared = point_radius * point_radius  # radius to point squared
        cos_phi0 = (point_radius_squared + radius_squared - range_squared) / (2 * radius * point_radius)
        sin_phi0 = np.sqrt(1 - cos_phi0 * cos_phi0)

        tmp1 = np.sqrt(4 * radius_squared * point_radius_squared - (
                point_radius_squared + radius_squared - range_squared) ** 2)
        phi0_rh = 2 * range / tmp1
        tmp2 = radar_sin_latitude * cos_phi0 + radar_cos_latitude * sin_phi0 * cos_azimuth
        tmp3 = 1 / np.sqrt(1 - tmp2 * tmp2)

        lat_ph = - tmp3 * radar_cos_latitude * sin_phi0 * sin_azimuth  # derivative of latitude with respect to azimuth
        lat_phi0 = tmp3 * (-radar_sin_latitude * sin_phi0 + radar_cos_latitude * cos_phi0 * cos_azimuth)
        lat_rh = lat_phi0 * phi0_rh  # derivative of latitude with respect to range

        sin_lat = radar_sin_latitude * cos_phi0 + radar_cos_latitude * sin_phi0 * cos_azimuth
        cos_lat_squared = 1 - sin_lat * sin_lat
        cos_lat = np.sqrt(cos_lat_squared)
        tmp4 = 1 / np.sqrt(1 - (sin_phi0 * sin_azimuth / cos_lat) ** 2)

        lon_ph = tmp4 * sin_phi0 * (
                    cos_azimuth * cos_lat + sin_azimuth * sin_lat * lat_ph) / cos_lat_squared  # derivative of longitude with respect to azimuth
        lon_rh = tmp4 * sin_azimuth * (
                cos_phi0 * phi0_rh * cos_lat + sin_phi0 * sin_lat * lat_rh) / cos_lat_squared  # derivative of longitude with respect to range

        ja = np.array([[lat_ph, lat_rh], [lon_ph, lon_rh]])  # Jacobian matrix

        # recalculation of the covariance matrix: ja * cov_radar * ja^T
        geodetic_covariance = ja @ np.array(cov_radar) @ ja.T

        return geodetic_covariance
