import numpy as np

from geodesy import EarthEllipsoid
import geometry_in_space
from CoordsCartesian3d import Cartesian3d


class StereoProjection:
    def __init__(self, center_geometry):
        self.point_geo_pos = np.zeros([3, ])
        self.point_stereo_pos = np.zeros([3, ])
        self.center_geo_pos = center_geometry.geo_pos
        self.center_geometry = center_geometry

    def set_geo_pos(self, point_geo_pos):
        self.point_geo_pos = point_geo_pos
        return

    def get_geo_pos(self):
        return self.point_geo_pos

    def set_stereo_pos(self, point_stereo_pos):
        self.point_stereo_pos = point_stereo_pos
        return

    def get_stereo_pos(self):
        return self.point_stereo_pos


    def geodetic_to_stereo(self, point_geo_pos):
        # Input: vector of geographic coordinates of the center center_geo_pos = [latitude, longitude, height]
        # latitude and longitude in radians, height in meters
        # vector of geographic coordinates point_geo_pos = [latitude, longitude, height]
        # latitude, longitude in radians, height in meters
        # Output: point coordinates in stereographic projection: point_stereo_pos = (x, y, height) in meters

        a2 = EarthEllipsoid.a2
        b2 = EarthEllipsoid.b2

        center_xyz_proj = self.center_geometry.xyz_proj
        xc, yc, zc = center_xyz_proj

        q0 = self.center_geometry.vector_to_the_nearest_pole
        n0 = self.center_geometry.normal

        # Cartesian coordinates of the point projection onto the ellipsoid
        point_geo_proj = [point_geo_pos[0], point_geo_pos[1], 0]
        cartesian3d_class = Cartesian3d()
        point_xyz_proj = np.array(cartesian3d_class.geodetic_to_cartesian3d(point_geo_proj))
        x0, y0, z0 = point_xyz_proj

        # equation of line (x + xc) / (x0 + xc) = (y + yc) / (y0 + yc) = (z + zc) / (z0 + zc) = t
        # equation of tangent plane x * xc + y * yc + z * zc * a^2 / b^2 = a^2
        t = 2 / (xc * x0 / a2 + yc * y0 / a2 + zc * z0 / b2 + 1)  # value of intersection parameter of line and plane

        # 3d Cartesian coordinates of the intersection point of a line and a plane, which is a stereographic point
        point_stereo_proj_3d = (point_xyz_proj + center_xyz_proj) * t - center_xyz_proj

        # vector from the center to the intersection point (parallel to the tangent plane)
        point_vector = point_stereo_proj_3d - center_xyz_proj

        # converting 3d coordinates of a point into 2d coordinates on the stereographic plane
        point_vector_length = np.linalg.norm(point_vector)
        if point_vector_length < self.center_geometry.zero_eps:  # we are too close to the projection center
            return [0, 0, point_geo_pos[2]]
        cos_alpha = np.dot(q0, point_vector) / point_vector_length

        sin_sign = np.sign(np.linalg.det([point_vector, q0, n0]))  # if the mixed product > 0, then the vectors form a right triple
        sin_alpha = sin_sign * np.sqrt(1 - cos_alpha * cos_alpha)

        point_stereo_pos = [point_vector_length * sin_alpha, point_vector_length * cos_alpha, point_geo_pos[2]]

        return point_stereo_pos


    def stereo_to_geodetic(self, point_stereo_pos):
        # Input: point coordinates in stereographic projection: point_stereo_pos = (x, y, height) in meters
        # vector of geographic coordinates of the center center_geo_pos
        # latitude and longitude in radians, height in meters
        # Output: vector of geographic coordinates point_geo_pos = (latitude, longitude, height)
        # latitude, longitude in radians, height in meters

        a2 = EarthEllipsoid.a2
        b2 = EarthEllipsoid.b2

        # Cartesian coordinates of the projection of the center of the system onto the ellipsoid
        center_xyz_proj = self.center_geometry.xyz_proj
        xc, yc, zc = center_xyz_proj

        q0 = self.center_geometry.vector_to_the_nearest_pole  # unit vector of north direction in the tangent plane
        n0 = self.center_geometry.normal  # unit normal vector to the ellipsoid at the center point

        # vector from center to intersection point
        x = point_stereo_pos[0]
        y = point_stereo_pos[1]

        point_vector_length = np.sqrt(x * x + y * y)
        alpha = np.atan2(x, y)

        # rotate the north direction vector in the tangent plane by an angle alpha,
        # we get a 3D vector that coincides in direction with the point direction vector in the tangent plane
        matrix = geometry_in_space.get_rotation_matrix(n0, -alpha)
        point_vector = point_vector_length * (matrix @ q0)

        # 3D projection coordinates
        point_stereo_proj_3d = point_vector + center_xyz_proj
        xs, ys, zs = point_stereo_proj_3d

        # equation of line (x + xc) / (xs + xc) = (y + yc) / (ys + yc) = (z + zc) / (zs + zc) = t
        # equation of ellipsoid x^2 / a^2 + y^2 / a^2 + z^2 / b^2 = 1
        # value of intersection parameter of line and ellipsoid
        denominator = ((xs + xc) ** 2 / a2 + (ys + yc) ** 2 / a2 + (zs + zc) ** 2 / b2)
        t = 2 * (xc * xs  / a2 + yc * ys / a2 + zc * zs / b2 + 1) / denominator

        point_xyz_pos = (point_stereo_proj_3d + center_xyz_proj) * t - np.array(center_xyz_proj)
        #point_xyz_pos = point_projection
        cartesian3d_class = Cartesian3d()
        point_geo_pos = cartesian3d_class.cartesian3d_to_geodetic(point_xyz_pos)
        point_geo_pos[2] = point_stereo_pos[2]

        return point_geo_pos


    def geodetic_to_stereo_covariance(self, point_geo_pos, geodetic_covariance):

        center_geometry = self.center_geometry

        cartesian3d_class = Cartesian3d()
        point_xyz_pos = cartesian3d_class.geodetic_to_cartesian3d(point_geo_pos)
        xyz_covariance = cartesian3d_class.geodetic_to_cartesian3d_covariance(point_geo_pos, geodetic_covariance)

        a2 = EarthEllipsoid.a2
        a2_b2 = EarthEllipsoid.a2_b2  # a2 / b2
        a4_b2 = EarthEllipsoid.a4_b2  # a4 / b2

        xc, yc, zc = center_geometry.xyz_pos  # Cartesian coordinates of the projection of the center of the projection onto the ellipsoid

        q0 = center_geometry.vector_to_the_nearest_pole  # unit vector of north direction in the tangent plane
        n0 = center_geometry.normal  # unit normal vector to the ellipsoid at the center point

        xa, ya, za = point_xyz_pos

        # derivatives of parameter t
        denominator = xc * xa + yc * ya + zc * za * a2_b2 + a2
        t = 2 * a2 / denominator

        denominator2 = denominator * denominator
        dt_dxa = -2 * a2 * xc / denominator2
        dt_dya = -2 * a2 * yc / denominator2
        dt_dza = -2 * a4_b2 * zc / denominator2

        # derivatives q = point_vector по xa, ya, za
        # point_vector = 2 * center_xyz_pos - (point_xyz_pos + center_xyz_pos) * t;
        dqx_dxa = t + (xa + xc) * dt_dxa
        dqx_dya = (xa + xc) * dt_dya
        dqx_dza = (xa + xc) * dt_dza

        dqy_dxa = (ya + yc) * dt_dxa
        dqy_dya = t + (ya + yc) * dt_dya
        dqy_dza = (ya + yc) * dt_dza

        dqz_dxa = (za + zc) * dt_dxa
        dqz_dya = (za + zc) * dt_dya
        dqz_dza = t + (za + zc) * dt_dza

        jacobian_matrix = np.array(
            [[dqx_dxa, dqx_dya, dqx_dza], [dqy_dxa, dqy_dya, dqy_dza], [dqz_dxa, dqz_dya, dqz_dza]])

        center_xyz_pos = np.array(center_geometry.xyz_pos)

        # s is a radius vector of the point_xyz_pos in 3d cartesian coordinates (3d array)
        s = (center_xyz_pos + point_xyz_pos) * t - center_xyz_pos
        s_norm = np.linalg.norm(s)

        # vector of partial derivatives of the function p
        d_p = np.array((q0 * s_norm ** 2 - s * np.dot(s, q0)) / s_norm ** 3)  # 3d vector

        # vector of partial derivatives of the alpha
        cos_alpha = np.dot(s, q0) / s_norm
        sin_alpha = np.sqrt(1 - cos_alpha * cos_alpha)
        if np.linalg.det([s, q0, n0]) < 0:
            sin_alpha = - sin_alpha
        d_alpha = -1 / sin_alpha * (jacobian_matrix.transpose() @ d_p)  # 3d vector

        # vector of partial derivatives of q_norm
        ds_norm = 1 / s_norm * (jacobian_matrix.transpose() @ s)  # 3d vector

        # derivatives of new coordinates with respect to old ones (required)
        dxS = sin_alpha * ds_norm + s_norm * cos_alpha * d_alpha  # 3d vector
        dyS = cos_alpha * ds_norm - s_norm * sin_alpha * d_alpha  # 3d vector

        # Jacobian matrix formation
        ja = np.array([dxS, dyS])  # 2x3 matrix

        # recalculate covariance matrix
        cov_stereo = ja @ xyz_covariance @ ja.T

        return cov_stereo
