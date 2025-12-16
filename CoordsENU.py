import numpy as np

import geometry_in_space
from CoordsCartesian3d import Cartesian3d

# ENU {East-North-Up) is a local tangent plane coordinate system used to describe positions relative
# to a reference point on the Earth's surface.

# ENU Coordinate Frame Definition
# At a given reference point (latitude \phi_c, longitude \lambda_c), the ENU coordinate system is defined as:
#    X-axis (East): Points eastward, tangent to the parallel.
#    Y-axis (North): Points northward, tangent to the meridian.
#    Z-axis (Up): Points away from the center of the Earth, perpendicular to the ellipsoid surface (i.e., local vertical).

class ENU:
    def __init__(self, center_geometry):
        self.center_geo_pos = center_geometry.geo_pos
        self.center_geometry = center_geometry


    def geodetic_to_enu(self, point_geo_pos):
        # Input: point_geo_pos = [latitude, longitude, height] are geodetic coordinates of the point
        # Output: point_enu_pos = [x, y, height] are ENU coordinates of the point

        center_xyz_pos = np.array(self.center_geometry.xyz_pos)

        point_xyz_pos = np.array(Cartesian3d.geodetic_to_cartesian3d(point_geo_pos))
        diff_vector = point_xyz_pos - center_xyz_pos
        rotation_matrix = geometry_in_space.get_rotation_matrix_for_enu(self.center_geometry.geo_pos)

        point_enu_pos = rotation_matrix @ diff_vector

        return point_enu_pos


    def enu_to_geodetic(self, point_enu_pos):
        # Input: point_enu_pos = [x, y, height] are ENU coordinates of the point
        # Output: point_geo_pos = [latitude, longitude, height] are geodetic coordinates of the point

        center_xyz_pos = np.array(self.center_geometry.xyz_pos)

        rotation_matrix = geometry_in_space.get_rotation_matrix_for_enu(self.center_geometry.geo_pos)
        diff_vector = rotation_matrix.T @ np.array(point_enu_pos)
        point_xyz_pos = diff_vector + center_xyz_pos

        point_geo_pos = np.array(Cartesian3d.cartesian3d_to_geodetic(point_xyz_pos))

        return point_geo_pos


    def geodetic_to_enu_covariance(self, point_geo_pos, geodetic_covariance):
        # conversion covariance from geodetic coordinates to scaledGeodetic

        xyz_covariance = Cartesian3d.geodetic_to_cartesian3d_covariance(point_geo_pos, geodetic_covariance)
        rotation_matrix = geometry_in_space.get_rotation_matrix_for_enu(self.center_geometry.geo_pos)

        # Jacobian matrix
        ja = rotation_matrix

        scaled_covariance = ja @ xyz_covariance @ ja.T

        return scaled_covariance
