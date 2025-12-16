import numpy as np

from geodesy import EarthEllipsoid
import geometry_in_space
from CoordScaledGeodetic import ScaledGeodeticCoordinates

class CartesianGeodesyCoordinates:
    def __init__(self, center_geometry):
        self.point_geo_pos = np.zeros([3, ])
        self.point_cart_local_pos = np.zeros([3, ])
        center_geo_pos = center_geometry.geo_pos
        center_geo_pos[1] = geometry_in_space.normalization_to_main_period(center_geo_pos[1])
        self.center_geo_pos = center_geo_pos
        self.center_geometry = center_geometry

    def set_geo_pos(self, point_geo_pos):
        self.point_geo_pos = point_geo_pos
        return

    def get_geo_pos(self):
        return self.point_geo_pos

    def set_polar_pos(self, point_cart_local_pos):
        self.point_cart_local_pos = point_cart_local_pos
        return

    def get_local_pos(self):
        return self.point_cart_local_pos


    def geodetic_to_local_cartesian(self, point_geo_pos):
        # Input: vector of geographic coordinates of the center center_geo_pos = [latitude, longitude, height]
        # vector of geographic coordinates point_geo_pos = [latitude, longitude, height]
        # latitude, longitude in radians, height in meters
        # Output: point coordinates in local Cartesian coordinates:
        # point_cart_local_pos = [parallel_length, meridian_length, height] in meters

        # bring the longitude of the system center to the interval [0, 2 * pi)
        # center_geo_pos[1] = geometry_on_plane.normalization_to_main_period(center_geo_pos[1])
        center_geo_pos = self.center_geo_pos

        # finding the y-coordinate from latitude
        latitude = point_geo_pos[0]
        meridian_length = EarthEllipsoid.meridian_to_meters_with_sign(center_geo_pos[0], latitude)
        # finding the x-coordinate
        longitude = geometry_in_space.normalization_to_main_period(point_geo_pos[1])  # bringing the longitude to [0, 2 * pi)

        # to cross the prime meridian
        d_longitude = longitude - center_geo_pos[1]
        if d_longitude > np.pi:
            d_longitude = d_longitude - 2 * np.pi
        if d_longitude < -np.pi:
            d_longitude = d_longitude + 2 * np.pi

        # convert longitude to local coordinate
        parallel_length = EarthEllipsoid.parallel_arc_length_with_sign(latitude, d_longitude)

        point_cart_local_pos = [parallel_length, meridian_length]

        return point_cart_local_pos


    def cartesian_local_to_geodetic(self, point_cart_local_pos):
        # Input: point coordinates in local Cartesian coordinates:
        # point_cart_local_pos = [parallel_length, meridian_length, height] in meters
        # vector of geographic coordinates of the center center_geo_pos = [latitude, longitude, height]
        # latitude and longitude in radians, height in meters
        # Output: vector of geographic coordinates point_geo_pos = [latitude, longitude, height]
        # latitude, longitude in radians [0, 2 * pi), height in meters

        center_geo_pos = self.center_geo_pos
        latitude_c = center_geo_pos[0]
        longitude_c = center_geo_pos[1]

        parallel_length = point_cart_local_pos[0]
        meridian_length = point_cart_local_pos[1]

        # calculating latitude
        d_latitude = EarthEllipsoid.meters_to_meridian_with_sign(latitude_c, meridian_length)
        latitude = latitude_c + d_latitude

        # crossing the pole
        if abs(latitude) > np.pi / 2:
            latitude = np.pi / 2 - latitude

        # calculating longitude
        d_longitude = EarthEllipsoid.meters_to_parallel(latitude, parallel_length)
        longitude = longitude_c + d_longitude
        longitude = geometry_in_space.normalization_to_main_period(longitude)

        point_geo_pos = [latitude, longitude]

        return point_geo_pos


    @staticmethod
    def geodetic_to_local_cartesian_covariance(geodetic_state):
        covariance = ScaledGeodeticCoordinates.geodetic_to_scaled_geodetic_covariance(geodetic_state)
        return covariance
