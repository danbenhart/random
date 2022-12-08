import matplotlib.pyplot as plt
import numpy as np
from scipy import optimize
from matplotlib.widgets import Slider

a = 1
b = 0.05
c = 7
start_angle = 0 * (2 * np.pi/360)
measured_degrees = 60 * (2 * np.pi/360)

full_radii = []
thetas = []

full_xs = []
full_ys = []


def fit_circle(x_list, y_list):

    x_m = np.mean(x_list)
    y_m = np.mean(y_list)

    def calc_R(xc, yc):
        """ calculate the distance of each 2D points from the center (xc, yc) """
        return np.sqrt((x_list - xc) ** 2 + (y_list - yc) ** 2)

    def f_2(c):
        """ calculate the algebraic distance between the data points and the mean circle centered at c=(xc, yc) """
        Ri = calc_R(*c)
        return Ri - Ri.mean()

    center_estimate = x_m, y_m
    # print(optimize.least_squares(f_2, center_estimate))
    center_2 = optimize.least_squares(f_2, center_estimate)['x']

    xc_2, yc_2 = center_2
    Ri_2 = calc_R(*center_2)
    R_2 = Ri_2.mean()
    residu_2 = sum((Ri_2 - R_2) ** 2)

    return xc_2, yc_2, R_2


# actual_xs = []
# actual_ys = []
#
# for theta in np.arange(0, 2 * np.pi, 2 * np.pi/360):
#     r = a + b * np.sin(c * theta)
#     full_radii.append(r)
#     thetas.append(theta)
#     x = r * np.cos(theta)
#     y = r * np.sin(theta)
#     full_xs.append(x)
#     full_ys.append(y)
#     if start_angle < theta < start_angle + measured_degrees:
#         actual_xs.append(x)
#         actual_ys.append(y)


def measured_stuff(nom_diameter=a,
                   nom_amplitude=b,
                   nom_wavelength=c,
                   start_angle=start_angle,
                   measured_degrees=measured_degrees):

    actual_xs = []
    actual_ys = []

    for theta in np.arange(0, 2 * np.pi, 2 * np.pi / 360):
        r = a + b * np.sin(c * theta)
        full_radii.append(r)
        thetas.append(theta)
        x = r * np.cos(theta)
        y = r * np.sin(theta)
        full_xs.append(x)
        full_ys.append(y)
        if start_angle < theta < start_angle + measured_degrees:
            actual_xs.append(x)
            actual_ys.append(y)

    measured_radii = []
    measured_thetas = []

    nominal_xs = []
    nominal_ys = []

    for degree in np.arange(0 + start_angle, measured_degrees + start_angle, 2 * np.pi/360):
        measured_radii.append(a)
        measured_thetas.append(degree)
        nominal_xs.append(a * np.cos(degree))
        nominal_ys.append(a * np.sin(degree))

    # fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
    # ax.plot(thetas, full_radii)
    # ax.plot(measured_thetas, measured_radii, color='red')
    # # ax.set_rmax(2)
    # ax.set_rticks([0.5, 1, 1.5, 2])  # Less radial ticks
    # ax.set_rlabel_position(-22.5)  # Move radial labels away from plotted line
    # ax.grid(True)
    #
    # ax.set_title("A line plot on a polar axis", va='bottom')
    # plt.show()

    actual_full_x, actual_full_y, actual_full_r = fit_circle(full_xs, full_ys)
    calc_x, calc_y, calc_r = fit_circle(actual_xs, actual_ys)
    # print(calc_r)

    nominal_full_circle = plt.Circle((actual_full_x, actual_full_y),
                                    actual_full_r,
                                    color='g',
                                    fill=False,
                                    label='nominal_circle')
    calc_circle = plt.Circle((calc_x, calc_y),
                             calc_r,
                             color='y',
                             fill=False,
                             label='calculated_circle')

    return actual_xs, actual_ys, nominal_xs, nominal_ys, nominal_full_circle, calc_circle


fig, ax = plt.subplots()
plt.subplots_adjust(left=0.5)
ax_start = plt.axes([0.1, 0.25, 0.0225, 0.63])
ax_dgrs_arc = plt.axes([0.2, 0.25, 0.0225, 0.63])
ax.set_aspect('equal', adjustable='box')

start_angle_slider = Slider(label='start_angle',
                            ax=ax_start,
                            valmin=0,
                            valmax=180,
                            valinit=0,
                            orientation='vertical')

dgrs_arc_slider = Slider(label='degrees of arc',
                         ax=ax_dgrs_arc,
                         valmin=0,
                         valmax=180,
                         valinit=60,
                         orientation='vertical')

actual_xs, actual_ys, nominal_xs, nominal_ys, nominal_full_circle, calc_circle = measured_stuff()


ax.plot(full_xs, full_ys, label='actual_circle', color='r')
ax.plot(nominal_xs, nominal_ys, label='nominal_radius', color='b')
ax.plot(actual_xs, actual_ys, label='measured_points', color='k')
ax.add_patch(nominal_full_circle)
ax.add_patch(calc_circle)
ax.legend(loc='upper left')


def update(val):
    ax.cla()
    start_angle = start_angle_slider.val * (2 * np.pi/360)
    degrees_arc = dgrs_arc_slider.val * (2 * np.pi/360)

    actual_xs, actual_ys, nominal_xs, nominal_ys, nominal_full_circle, calc_circle = measured_stuff(start_angle=start_angle,
                                                                                                    measured_degrees=degrees_arc)

    ax.plot(full_xs, full_ys, label='actual_circle', color='r')
    ax.plot(nominal_xs, nominal_ys, label='nominal_radius', color='b')
    ax.plot(actual_xs, actual_ys, label='measured_points', color='k')
    ax.add_patch(nominal_full_circle)
    ax.add_patch(calc_circle)
    ax.legend(loc='upper left')
    fig.canvas.draw()


start_angle_slider.on_changed(update)
dgrs_arc_slider.on_changed(update)

plt.show()

