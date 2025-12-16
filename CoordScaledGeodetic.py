import numpy as np

from geodesy import EarthEllipsoid


class ScaledGeodeticCoordinates:
    def __init__(self, center_geo_pos = None):
        self.point_geo_pos = np.zeros([3, ])
        self.point_scaled_geo_pos = np.zeros([3, ])
        self.center_geo_pos = center_geo_pos

    def set_geo_pos(self, point_geo_pos):
        self.point_geo_pos = point_geo_pos
        return

    def get_geo_pos(self):
        return self.point_geo_pos

    def set_polar_pos(self, point_scaled_geo_pos):
        self.point_scaled_geo_pos = point_scaled_geo_pos
        return

    def get_local_pos(self):
        return self.point_scaled_geo_pos


    def geodetic_to_scaled_geodetic(self, point_geo_pos):
        # conversion from geographic coordinates to scaledGeodetic and scaledGeodetic0 coordinates

        radius = EarthEllipsoid.average_radius

        # by latitude
        latitude = point_geo_pos[0]
        y = latitude * radius  # distance from the equator along the meridian on the sphere (in meters)

        # by longitude
        if self.center_geo_pos is None:
            parallel_lat = latitude  # recalculation in scaledGeodetic0
        else:
            parallel_lat = self.center_geo_pos[0]  # recalculation in scaledGeodetic

        longitude = point_geo_pos[1]
        scale_x = radius * np.cos(parallel_lat)  # radius of parallel
        x = longitude * scale_x  # distance from Greenwich along the arc of parallel on the sphere (in meters)

        point_scaled_geo_pos = [x, y, point_geo_pos[2]]

        return point_scaled_geo_pos


    def scaled_geodetic_to_geodetic(self, point_scaled_geo_pos):
        # conversion from scaledGeodetic and scaledGeodetic0 coordinates to geodetic coordinates

        radius = EarthEllipsoid.average_radius

        # by latitude
        y = point_scaled_geo_pos[1]  # arc length from the equator along the meridian (in meters)
        latitude = y / radius

        if self.center_geo_pos is None:
            parallel_lat = latitude   # recalculation from scaledGeodetic0
        else:
            parallel_lat = self.center_geo_pos[0]  # recalculation from scaledGeodetic

        # by longitude
        scale_x = radius * np.cos(parallel_lat)  # radius of parallel
        x = point_scaled_geo_pos[0]  # distance from Greenwich along the arc of parallel (in meters)
        longitude = x / scale_x

        point_geo_pos = [latitude, longitude, point_scaled_geo_pos[2]]

        return point_geo_pos


    @staticmethod
    def geodetic_to_scaled_geodetic_covariance(geodetic_state):
        # conversion covariance from geodetic coordinates to scaledGeodetic

        point_geo_pos = geodetic_state.X
        geodetic_covariance = geodetic_state.P
        radius = EarthEllipsoid.average_radius
        latitude = point_geo_pos[0]

        scale_x = radius * np.cos(latitude)  # radius of parallel
        scale_y = radius  # mean radius of the earth (and meridian)

        # Jacobian matrix
        ja = np.array([[0, scale_x],
                       [scale_y, 0]])

        scaled_covariance = ja @ geodetic_covariance @ ja.T

        return scaled_covariance
