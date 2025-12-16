import matplotlib.pyplot as plt
import numpy as np


def plot_errors(error_array, error_name):
    # error chart

    fig, ax = plt.subplots()
    ax.plot(error_array)

    ax.set(xlabel='point number', ylabel='error value',
           title=error_name)
    ax.grid()
    plt.show()

    return


def plot_coord_errors_for_radar(graph_on, point_geo_pos, err_geo, err_polar, point_range):

    number_of_points = len(point_geo_pos)

    if graph_on:
        fig1, plt1 = plt.subplots(3, 1)
        point_numbers = np.array(range(number_of_points))
        names_geo = ["latitude", "longitude", "altitude"]
        for k in range(3):
            point_geo_pos = np.array(point_geo_pos)
            plt1[k].plot(point_geo_pos[:, k], err_geo[:, k], label=names_geo[k])
            plt1[k].legend()

            plt1[k].set_xlabel("point number")
            plt1[k].set_ylabel("error")
            plt1[k].grid()

        plt.subplots_adjust(hspace=0.5)
        fig2, plt2 = plt.subplots(4, 1)
        names_polar = ["azimuth", "range", "height"]
        for k in range(3):
            plt2[k].plot(point_numbers, err_polar[:, k], label=names_polar[k])
            plt2[k].legend()

            plt2[k].set_xlabel("point number")
            plt2[k].set_ylabel("error")
            plt2[k].grid()

        plt2[3].plot(point_numbers, point_range, label=range)
        plt2[3].legend()

        plt2[3].set_xlabel("point number")
        plt2[3].set_ylabel("error")
        plt2[3].grid()

        plt.subplots_adjust(hspace=0.5)

        plt.show()

    return


def plot_coord_cart3d_errors(graph_on, point_geo_pos, err_geo, err_cart3d):

    number_of_points = len(point_geo_pos)

    if graph_on:
        fig1, plt1 = plt.subplots(3, 1)
        point_numbers = np.array(range(number_of_points))
        names_geo_err = ["latitude error", "longitude error", "altitude error"]
        for k in range(3):
           point_geo_pos = np.array(point_geo_pos)
           plt1[k].plot(point_geo_pos[:, k], err_geo[:, k], label = names_geo_err[k])
           plt1[k].legend()

           plt1[k].grid()

        plt.subplots_adjust(hspace=0.5)

        fig2, plt2 = plt.subplots(3, 1)
        names_xyz_err = ["x error", "y error", "z error"]
        for k in range(3):
            plt2[k].plot(point_numbers, err_cart3d[:, k], label=names_xyz_err[k])
            plt2[k].legend()

            plt2[k].grid()

        plt.subplots_adjust(hspace=0.5)

        plt.show()

    return


def plot_coord_stereo_errors(graph_on, stereo_range, err_geo, err_stereo):

    number_of_points = len(stereo_range)

    if graph_on:
        fig1, plt1 = plt.subplots(3, 1)
        point_numbers = np.array(range(number_of_points))
        names_geo_err = ["latitude error", "longitude error"]

        plt1[0].plot(point_numbers, err_geo[:, 0], label = names_geo_err[0])
        plt1[1].plot(point_numbers, err_geo[:, 1], label=names_geo_err[1])
        plt1[2].plot(point_numbers, stereo_range, label="dist from center")
        for k in range(3):
            plt1[k].legend()
            plt1[k].grid()

        plt.subplots_adjust(hspace=0.5)
        fig2, plt2 = plt.subplots(3, 1)
        names_xyz_err = ["x error", "y error"]

        plt2[0].plot(point_numbers, err_stereo[:, 0], label=names_xyz_err[0])
        plt2[1].plot(point_numbers, err_stereo[:, 1], label=names_xyz_err[1])
        plt2[2].plot(point_numbers, stereo_range, label="dist from center")
        for k in range(3):
            plt2[k].legend()
            plt2[k].grid()

        plt.subplots_adjust(hspace=0.5)
        plt.show()

    return


def plot_coord_radar_errors(graph_on, slant_range, err_geo, err_polar):

    number_of_points = len(slant_range)

    if graph_on:
        fig1, plt1 = plt.subplots(3, 1)
        point_numbers = np.array(range(number_of_points))
        names_geo_err = ["latitude error", "longitude error"]

        plt1[0].plot(point_numbers, err_geo[:, 0], label = names_geo_err[0])
        plt1[1].plot(point_numbers, err_geo[:, 1], label = names_geo_err[1])
        plt1[2].plot(point_numbers, slant_range, label="dist from center")
        for k in range(3):
            plt1[k].legend()
            plt1[k].grid()

        plt.subplots_adjust(hspace=0.5)

        fig2, plt2 = plt.subplots(3, 1)
        names_xyz_err = ["azimuth error", "range error"]

        plt2[0].plot(point_numbers, err_polar[:, 0], label=names_xyz_err[0])
        plt2[1].plot(point_numbers, err_polar[:, 1], label=names_xyz_err[1])
        plt2[2].plot(point_numbers, slant_range, label="dist from center")
        for k in range(3):
            plt2[k].legend()
            plt2[k].grid()

        plt.subplots_adjust(hspace=0.5)


        plt.show()

    return


def plot_coord_scaled_geo_errors(graph_on, point_geo_pos, err_geo, err_local):

    number_of_points = len(point_geo_pos)

    if graph_on:
        fig1, plt1 = plt.subplots(2, 1)
        point_numbers = np.array(range(number_of_points))
        names_geo_err = ["latitude error", "longitude error"]
        for k in range(2):
           point_geo_pos = np.array(point_geo_pos)
           plt1[k].plot(point_geo_pos[:, k], err_geo[:, k], label = names_geo_err[k])
           plt1[k].legend()

           plt1[k].grid()

        plt.subplots_adjust(hspace=0.5)

        fig2, plt2 = plt.subplots(2, 1)
        names_xyz_err = ["x error", "y error"]
        for k in range(2):
            plt2[k].plot(point_numbers, err_local[:, k], label=names_xyz_err[k])
            plt2[k].legend()

            plt2[k].grid()

        plt.subplots_adjust(hspace=0.5)

        plt.show()

    return


def plot_coord_cart_local_errors(graph_on, dist_from_center, err_geo, err_local):

    number_of_points = len(dist_from_center)

    if graph_on:
        fig1, plt1 = plt.subplots(3, 1)
        point_numbers = np.array(range(number_of_points))
        names_geo_err = ["latitude error", "longitude error"]

        plt1[0].plot(point_numbers, err_geo[:, 0], label = names_geo_err[0])
        plt1[1].plot(point_numbers, err_geo[:, 1], label=names_geo_err[1])
        plt1[2].plot(point_numbers, dist_from_center, label="dist from center")
        for k in range(3):
            plt1[k].legend()
            plt1[k].grid()

        plt.subplots_adjust(hspace=0.5)
        fig2, plt2 = plt.subplots(3, 1)
        names_xyz_err = ["x error", "y error"]

        plt2[0].plot(point_numbers, err_local[:, 0], label=names_xyz_err[0])
        plt2[1].plot(point_numbers, err_local[:, 1], label=names_xyz_err[1])
        plt2[2].plot(point_numbers, dist_from_center, label="dist from center")
        for k in range(3):
            plt2[k].legend()
            plt2[k].grid()

        plt.subplots_adjust(hspace=0.5)
        plt.show()

    return

##########################################
def plot_coord_enu_errors(graph_on, enu_range, err_geo, err_enu):

    number_of_points = len(enu_range)

    if graph_on:
        fig1, plt1 = plt.subplots(3, 1)
        point_numbers = np.array(range(number_of_points))
        names_geo_err = ["latitude error", "longitude error"]

        plt1[0].plot(point_numbers, err_geo[:, 0], label = names_geo_err[0])
        plt1[1].plot(point_numbers, err_geo[:, 1], label=names_geo_err[1])
        plt1[2].plot(point_numbers, enu_range, label="dist from center")
        for k in range(3):
            plt1[k].legend()
            plt1[k].grid()

        plt.subplots_adjust(hspace=0.5)
        fig2, plt2 = plt.subplots(3, 1)
        names_xyz_err = ["x error", "y error"]

        plt2[0].plot(point_numbers, err_enu[:, 0], label=names_xyz_err[0])
        plt2[1].plot(point_numbers, err_enu[:, 1], label=names_xyz_err[1])
        plt2[2].plot(point_numbers, enu_range, label="dist from center")
        for k in range(3):
            plt2[k].legend()
            plt2[k].grid()

        plt.subplots_adjust(hspace=0.5)
        plt.show()

    return
