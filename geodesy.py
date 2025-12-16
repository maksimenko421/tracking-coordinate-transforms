import numpy as np

import general_settings

def calc_parameters(a, e2):
    # parameters for calculating the radius of curvature $N$ of the meridian
    # the coefficients of the Taylor series expansion of the M = a * (1 - e^2) / sqrt((1 - e^2 * sin(latitude)^2)^3)
    m0 = a * (1 - e2)
    m2 = 1.5 * e2 * m0
    m4 = 1.25 * e2 * m2
    m6 = 7 / 6 * e2 * m4
    m8 = 9 / 8 * e2 * m6
    m10 = 1.1 * e2 * m8
    m_parameters = [m0, m2, m4, m6, m8, m10]

    # parameters for calculating the radius of curvature of a parallel
    # the coefficients of the Taylor series expansion of the N = a / sqrt(1 - e^2 * sin(latitude)^2)
    n0 = a
    n2 = 0.5 * e2 * n0
    n4 = 0.75 * e2 * n2
    n6 = 5 / 6 * e2 * n4
    n8 = 7 / 8 * e2 * n6
    n10 = 0.9 * e2 * n8
    n_parameters = [n0, n2, n4, n6, n8, n10]

    # parameters for calculating the meridian arc length in meters
    a0 = m0 + m2 / 2 + 3 / 8 * m4 + 5 / 16 * m6 + 35 / 128 * m8
    a2 = m2 / 4 + m4 / 4 + 15 / 64 * m6 + 7 / 32 * m8
    a4 = m4 / 32 + 3 / 64 * m6 + 7 / 128 * m8
    a6 = m6 / 192 + m8 / 96
    a_parameters = [a0, a2, a4, a6]

    # parameters for calculating latitude by the length of the meridian arc
    b0 = a0
    b2 = a2 / a0 * (1 + a4 / a0 - a2 * a2 / (2 * a0 * a0))
    b4 = a2 * a2 / (a0 * a0) - a4 / a0
    b6 = a6 / a0 - 3 * a2 / a0 * (a4 / a0 - a2 * a2 / (2 * a0 * a0))
    b_parameters = [b0, b2, b4, b6]

    return m_parameters, n_parameters, a_parameters, b_parameters


class EarthEllipsoid:
    a = general_settings.Ellipsoid["a"]
    a2 = a * a
    a4 = a2 * a2
    b = general_settings.Ellipsoid["b"]
    b2 = b * b
    b4 = b2 * b2
    c = np.sqrt(a2 - b2)  # focal length
    e = c / a  # eccentricity
    e2 = 1 - b2 / a2 # eccentricity squared
    f = (a - b) / a  # flattening
    alpha = 1 / f  # inverse flattening
    average_radius = 6371000
    a2_b2 = a2 / b2  # a2 / b2
    a4_b2 = a4 / b2  # a4 / b2
    a4_b4 = a2_b2 * a2_b2   # (a2 / b2) ** 2

    m_parameters, n_parameters, a_parameters, b_parameters = calc_parameters(a, e2)


    @staticmethod
    def meridian_curvature_radius(latitude):
        # radius of curvature of the meridian
        # M = a * (1 - e^2) / sqrt((1 - e^2 * sin(latitude)^2)^3)
        # Morozov "Course of spherical geodesy" p.20

        a = EarthEllipsoid.a
        e2 = EarthEllipsoid.e2
        sin_lat = np.sin(latitude)

        meridian_radius = a * (1 - e2) / np.sqrt((1 - e2 * sin_lat * sin_lat) ** 3)

        return meridian_radius


    @staticmethod
    def radius_of_curvature_of_the_prime_vertical(latitude):
        # radius of curvature of the prime vertical
        # N = a / sqrt(1 - e^2 * sin(latitude)^2)
        # Morozov "Course of spherical geodesy" p.20

        a = EarthEllipsoid.a
        e2 = EarthEllipsoid.e2
        sin_lat = np.sin(latitude)

        vertical_radius = a / np.sqrt(1 - e2 * sin_lat * sin_lat)

        return vertical_radius


    @staticmethod
    def meridian_to_meters_from_equator(latitude):
        # calculation of the length of the meridian arc from zero latitude to latitude
        # Morozov "Course of spherical geodesy" p.28
        # maximum error is 0.0001 meters

        if latitude < 0:
            latitude = -latitude

        [a0, a2, a4, a6] = EarthEllipsoid.a_parameters  # parameters are already divided by constants 1, 2, 4, 6
        meridian_length = a0 * latitude - a2 * np.sin(2 * latitude) + a4 * np.sin(4 * latitude) - a6 * np.sin(6 * latitude)

        return meridian_length


    @staticmethod
    def meters_to_meridian_from_equator(meridian_length):
        # calculation of latitude (positive) by the length of the meridian arc from the equator
        # Morozov "Course of Spherical Geodesy" p.29
        # the maximum error is 0.005 seconds

        [b0, b2, b4, b6] = EarthEllipsoid.b_parameters

        beta = meridian_length / b0
        latitude = beta + b2 * np.sin(2 * beta) + b4 * np.sin(4 * beta) + b6 * np.sin(6 * beta)

        return latitude


    @staticmethod
    def meridian_to_meters_with_sign(latitude1, latitude2):
        # calculation of the meridian arc length from latitude latitude1 to latitude latitude2
        # Morozov "Course of Spherical Geodesy" p.31
        # maximum error is 0.002 meters for small arc differences (latitude differences up to 30')

        meridian_length1 = EarthEllipsoid.meridian_to_meters_from_equator(latitude1)
        if latitude1 < 0:
            meridian_length1 = - meridian_length1
        meridian_length2 = EarthEllipsoid.meridian_to_meters_from_equator(latitude2)
        if latitude2 < 0:
            meridian_length2 = - meridian_length2

        meridian_length = meridian_length2 - meridian_length1  # the latitudes are in one semi-sphere

        return meridian_length


    @staticmethod
    def meters_to_meridian_with_sign(latitude, meridian_length):
        # calculation of the difference in latitudes by the length of the meridian arc from the latitude
        # output lat_diff is always positive
        # Morozov "Course of spherical geodesy" p. 31
        # the maximum error is 0.002 meters for small differences in arcs (differences in latitude up to 30')

        sign = 1
        if latitude < 0:
            sign = -1
            latitude = -latitude

        meridian_latitude = EarthEllipsoid.meridian_to_meters_from_equator(latitude)
        latitude_common = EarthEllipsoid.meters_to_meridian_from_equator(meridian_length + meridian_latitude)
        lat_diff = latitude_common - latitude

        return sign * lat_diff


    @staticmethod
    def meters_to_parallel(latitude, parallel_length):
        # calculation of the difference longitude_diff by the length of the arc of the parallel parallel_length
        # Morozov "Course of spherical geodesy" p. 32

        if latitude < 0:
            latitude = -latitude

        vertical_radius = EarthEllipsoid.radius_of_curvature_of_the_prime_vertical(latitude)
        longitude_diff = parallel_length / (vertical_radius * np.cos(latitude))

        return longitude_diff


    @staticmethod
    def parallel_arc_length_with_sign(latitude, longitude_diff):
        # calculation of the length of an arc of parallel at a given latitude
        # Morozov "Course of Spherical Geodesy" p.32

        vertical_radius = EarthEllipsoid.radius_of_curvature_of_the_prime_vertical(latitude)
        parallel_arc_length = vertical_radius * np.cos(latitude) * longitude_diff

        return parallel_arc_length


    @staticmethod
    def geocentric_radius(latitude):
        # in meters
        a2 = EarthEllipsoid.a2
        b2 = EarthEllipsoid.b2

        cos_lat = np.cos(latitude)
        sin_lat = np.sin(latitude)

        numerator = (a2 * cos_lat) ** 2 + (b2 * sin_lat) ** 2
        denominator = a2 * cos_lat ** 2 + b2 * sin_lat ** 2

        geocentric_radius = np.sqrt(numerator / denominator)

        return geocentric_radius


    @staticmethod
    def geodetic_string_to_radians(str_latitude, str_longitude):
        # convert the latitude and longitude specified by the string to numbers
        # input latitude, longitude - strings
        # output longitude in [0, 2*pi]  (in radians)
        # output latitude in [-pi/2, pi/2]  (in radians)

        latitude = int(str_latitude)
        longitude = int(str_longitude)

        # finding longitude
        sign_longitude = np.sign(longitude)
        longitude = np.abs(longitude)

        seconds = longitude % 100
        minutes = round((longitude - seconds) / 100) % 100
        degrees = round((longitude - seconds - 100 * minutes) / 10000)

        longitude = (degrees + minutes / 60 + seconds / 60 / 60) * np.pi / 180

        if sign_longitude < 0:
            longitude = 2 * np.pi - longitude

        # finding latitude
        sign_latitude = np.sign(latitude)
        latitude = np.abs(latitude)

        seconds = latitude % 100
        minutes = round((latitude - seconds) / 100) % 100
        degrees = round((latitude - seconds - 100 * minutes) / 10000)

        latitude = sign_latitude * (degrees + minutes / 60 + seconds / 60 / 60) * np.pi / 180

        return latitude, longitude
