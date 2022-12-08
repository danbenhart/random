import matplotlib.pyplot as plt
import numpy as np
from scipy import optimize
from matplotlib.widgets import Slider

a = 1
b = 0.05
c = 3

b2 = 0.05
c2 = 5

start_angle = 0 * (2 * np.pi / 360)
measured_degrees = 60 * (2 * np.pi / 360)

# full_radii = []
# thetas = []
#
# full_xs = []
# full_ys = []


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

    # center_estimate = x_m, y_m
    center_estimate = 0, 0
    # print(optimize.least_squares(f_2, center_estimate))
    center_2 = optimize.least_squares(f_2, center_estimate)['x']

    xc_2, yc_2 = center_2
    Ri_2 = calc_R(*center_2)
    R_2 = Ri_2.mean()
    residu_2 = sum((Ri_2 - R_2) ** 2)

    return xc_2, yc_2, R_2


def measured_stuff(nom_diameter=a,
                   form_1_amplitude=b,
                   form_2_amplitude=b2,
                   form_1_wavelength=c,
                   form_2_wavelength=c2,
                   start_angle=start_angle,
                   measured_degrees=measured_degrees):
    actual_xs = []
    actual_ys = []

    full_radii = []
    thetas = []

    full_xs = []
    full_ys = []

    for theta in np.arange(0, 2 * np.pi, 2 * np.pi / 360):
        r = nom_diameter + form_1_amplitude * np.sin(form_1_wavelength * theta)
        r = r * (nom_diameter + form_2_amplitude * np.sin(form_2_wavelength * theta))
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

    for degree in np.arange(0 + start_angle, measured_degrees + start_angle, 2 * np.pi / 360):
        measured_radii.append(a)
        measured_thetas.append(degree)
        nominal_xs.append(a * np.cos(degree))
        nominal_ys.append(a * np.sin(degree))

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

    return full_xs, full_ys, actual_xs, actual_ys, nominal_xs, nominal_ys, nominal_full_circle, calc_circle


fig, ax = plt.subplots()
plt.subplots_adjust(left=0.5)

ax_start = plt.axes([0.05, 0.25, 0.0225, 0.63])
ax_dgrs_arc = plt.axes([0.1, 0.25, 0.0225, 0.63])
ax_form_1_wavelength = plt.axes([0.15, 0.25, 0.0225, 0.63])
ax_form_1_amp = plt.axes([0.2, 0.25, 0.0225, 0.63])
ax_form_2_wavelength = plt.axes([0.25, 0.25, 0.0225, 0.63])
ax_form_2_amp = plt.axes([0.3, 0.25, 0.0225, 0.63])

ax.set_aspect('equal', adjustable='box')

start_angle_slider = Slider(label='start_angle',
                            ax=ax_start,
                            valmin=0,
                            valmax=180,
                            valinit=0,
                            valstep=1,
                            orientation='vertical')

dgrs_arc_slider = Slider(label='degrees of arc',
                         ax=ax_dgrs_arc,
                         valmin=0,
                         valmax=180,
                         valinit=60,
                         valstep=1,
                         orientation='vertical')

form_1_wavelength_slider = Slider(label='form_1\nwavelengths',
                                  ax=ax_form_1_wavelength,
                                  valmin=1,
                                  valmax=10,
                                  valinit=c,
                                  valstep=1,
                                  orientation='vertical')

form_1_amplitude_slider = Slider(label='form_1\namplitude',
                                 ax=ax_form_1_amp,
                                 valmin=0,
                                 valmax=0.1,
                                 valinit=0.05,
                                 # valstep=1,
                                 orientation='vertical')

form_2_wavelength_slider = Slider(label='form_2\nwavelengths',
                                  ax=ax_form_2_wavelength,
                                  valmin=1,
                                  valmax=180,
                                  valinit=c,
                                  valstep=1,
                                  orientation='vertical')

form_2_amplitude_slider = Slider(label='form_2\namplitude',
                                 ax=ax_form_2_amp,
                                 valmin=0,
                                 valmax=0.1,
                                 valinit=0.05,
                                 # valstep=1,
                                 orientation='vertical')


full_xs, full_ys, actual_xs, actual_ys, nominal_xs, nominal_ys, nominal_full_circle, calc_circle = measured_stuff()

ax.plot(full_xs, full_ys, label='actual_circle', color='r')
ax.plot(nominal_xs, nominal_ys, label='nominal_radius', color='b')
ax.plot(actual_xs, actual_ys, label='measured_points', color='k')
ax.add_patch(nominal_full_circle)
ax.add_patch(calc_circle)
ax.legend(loc='upper left')


def update(val):
    ax.cla()
    start_angle = start_angle_slider.val * (2 * np.pi / 360)
    degrees_arc = dgrs_arc_slider.val * (2 * np.pi / 360)
    form_1_wavelengths = int(form_1_wavelength_slider.val)
    form_1_amp = form_1_amplitude_slider.val
    form_2_wavelengths = int(form_2_wavelength_slider.val)
    form_2_amp = form_2_amplitude_slider.val

    full_xs, full_ys, actual_xs, actual_ys, nominal_xs, nominal_ys, nominal_full_circle, calc_circle = measured_stuff(
        start_angle=start_angle,
        measured_degrees=degrees_arc,
        form_1_amplitude=form_1_amp,
        form_1_wavelength=form_1_wavelengths,
        form_2_amplitude=form_2_amp,
        form_2_wavelength=form_2_wavelengths)

    ax.plot(full_xs, full_ys, label='actual_circle', color='r')
    ax.plot(nominal_xs, nominal_ys, label='nominal_radius', color='b')
    ax.plot(actual_xs, actual_ys, label='measured_points', color='k')
    ax.add_patch(nominal_full_circle)
    ax.add_patch(calc_circle)
    ax.legend(loc='upper left')
    fig.canvas.draw()


start_angle_slider.on_changed(update)
dgrs_arc_slider.on_changed(update)
form_1_wavelength_slider.on_changed(update)
form_1_amplitude_slider.on_changed(update)
form_2_wavelength_slider.on_changed(update)
form_2_amplitude_slider.on_changed(update)

plt.show()

