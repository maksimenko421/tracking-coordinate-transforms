import numpy as np

from geodesy import EarthEllipsoid
import geometry_in_space

# Earth - Centered, Earth - Fixed(ECEF)
# Cartesian 3D coordinates: [x, y, z] in [meters, meters, meters]:
# The coordinate center is at the center of the Earth, the oz axis passes through the North Pole,
# the ox axis is in the plane of the equator and passed through 0th meridian, the oy axis is in the equator plane
# and is perpendicular ox so, the ox. oy, oz are the right vector triple
#
# Geodetic coordinates [latitude, longitude, altitude] in [radians, radians, meters]
# the latitude is an acute angle formed by the normal to the surface of an ellipsoid and the plane of the equator,
# belongs to [0, pi/2] in the northern hemisphere and to [-pi/2,0] in the southern hemisphere.
# the longitude is a dihedral angle formed by the plane of the prime meridian and a given meridian,
# belongs to [0, pi] in the eastern hemisphere and to [-pi,0] in the western hemisphere.
# the altitude is the distance from the point to the surface of the ellipsoid

class Cartesian3d:
    def __init__(self):
        self.point_geo_pos = np.zeros([3, ])
        self.point_xyz_pos = np.zeros([3, ])

    def set_geo_pos(self, point_geo_pos):
        self.point_geo_pos = point_geo_pos
        return

    def get_geo_pos(self):
        return self.point_geo_pos

    def set_xyz_pos(self, point_xyz_pos):
        self.point_xyz_pos = point_xyz_pos
        return

    def get_xyz_pos(self):
        return self.point_xyz_pos


    @staticmethod
    def geodetic_to_cartesian3d(point_geo_pos):
        # Input: point_geo_pos = [latitude, longitude, height] are Geodetic coordinates of the point
        # Output: point_xyz_pos = [x, y, z] are Cartesian 3D coordinates of the point

        e2 = EarthEllipsoid.e2  # square of eccentricity
        latitude, longitude, height = point_geo_pos

        # radius of curvature of the prime vertical
        vertical_radius = EarthEllipsoid.radius_of_curvature_of_the_prime_vertical(latitude)

        x = (vertical_radius + height) * np.cos(latitude) * np.cos(longitude)
        y = (vertical_radius + height) * np.cos(latitude) * np.sin(longitude)
        z = ((1 - e2) * vertical_radius + height) * np.sin(latitude)

        point_xyz_pos = [x, y, z]

        return point_xyz_pos


    @staticmethod
    def cartesian3d_to_geodetic(point_xyz_pos):
        # Input: point_xyz_pos = [xa, ya, za] are Cartesian 3D coordinates of the point
        # Output: point_geo_pos = [latitude, longitude, altitude] are Geodetic coordinates of the point

        e2 = EarthEllipsoid.e2  # eccentricity squared
        x, y, z = point_xyz_pos

        # Find projection of the point onto ellipsoid
        x0, y0, z0 = geometry_in_space.project_onto_ellipsoid(point_xyz_pos)

        # Find needed values
        height = np.sqrt((x - x0) ** 2 + (y - y0) ** 2 + (z - z0) ** 2)
        longitude = np.atan2(y, x)
        longitude = geometry_in_space.normalization_to_main_period(longitude)
        latitude = np.atan((z0 / np.sqrt(x0 * x0 + y0 * y0)) / (1 - e2))

        point_geo_pos = [latitude, longitude, height]

        return point_geo_pos


    @staticmethod
    def geodetic_to_cartesian3d_covariance(point_geo_pos, geodetic_covariance, sigma2_h=1):
        a = EarthEllipsoid.a
        e2 = EarthEllipsoid.e2

        latitude, longitude, height = point_geo_pos
        sin_lat = np.sin(latitude)
        cos_lat = np.cos(latitude)
        sin_lon = np.sin(longitude)
        cos_lon = np.cos(longitude)

        # radius of curvature of the first vertical
        vertical_radius = EarthEllipsoid.radius_of_curvature_of_the_prime_vertical(latitude)

        general_height = vertical_radius + height

        # we find the derivatives and the Jacobian matrix
        dN_dlat = a * e2 * sin_lat * cos_lat / (1 - e2 * sin_lat * sin_lat) ** (3 / 2)

        dx_dlat = (dN_dlat * cos_lat - general_height * sin_lat) * cos_lon
        dx_dlon = -general_height * cos_lat * sin_lon
        dx_dhei = cos_lat * cos_lon

        dy_dlat = (dN_dlat * cos_lat - general_height * sin_lat) * sin_lon
        dy_dlon = general_height * cos_lat * cos_lon
        dy_dhei = cos_lat * sin_lon

        dz_dlat = dN_dlat * (1 - e2) * sin_lat + ((1 - e2) * vertical_radius + height) * cos_lat
        dz_dlon = 0
        dz_dhei = sin_lat

        # Jacobian matrix
        ja = np.array([[dx_dlat, dx_dlon, dx_dhei],
                       [dy_dlat, dy_dlon, dy_dhei],
                       [dz_dlat, dz_dlon, dz_dhei]])

        # recalculation of the covariance matrix
        cov_geod = np.zeros([3, 3])
        cov_geod[0:2, 0:2] = geodetic_covariance
        cov_geod[2, 2] = sigma2_h
        covariance_xyz = ja @ cov_geod @ ja.T

        return covariance_xyz


class PointGeometry:
    def __init__(self, point_geo_pos, needed_type="extended"):
        [latitude, longitude, height] = point_geo_pos
        self.geo_pos = np.array(point_geo_pos)
        cartesian3d_class = Cartesian3d()
        self.xyz_pos = np.array(cartesian3d_class.geodetic_to_cartesian3d(self.geo_pos))
        self.geo_proj = np.array([latitude, longitude, 0])
        self.xyz_proj = np.array(cartesian3d_class.geodetic_to_cartesian3d(self.geo_proj))
        if needed_type != "simple":
            self.normal = geometry_in_space.normal_to_ellipsoid(self.xyz_proj)
            self.vector_to_the_nearest_pole = geometry_in_space.vector_towards_the_nearest_pole(self.xyz_proj)
            self.normal_vector_to_meridian_plane = geometry_in_space.normal_to_the_meridian_plane(self.xyz_proj)
            self.radius = np.linalg.norm(self.xyz_pos)  # radius of the sphere at the radar point
            self.geocentric_radius = EarthEllipsoid.geocentric_radius(latitude)
            self.sin_latitude = np.sin(latitude)
            self.cos_latitude = np.cos(latitude)
            self.zero_eps = 1e-7
