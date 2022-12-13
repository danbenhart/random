import matplotlib.pyplot as plt
import numpy as np
from scipy import optimize
from matplotlib.widgets import Slider

init_nom_radius = 1
init_form_1_amp = 0.05
init_form_1_wave = 3

init_form_2_amp = 0.05
init_form_2_wave = 5

init_start_angle = 0 * (2 * np.pi / 360)
init_measured_degrees = 60 * (2 * np.pi / 360)

form_1_max_amp = .001     # worst case scenario for lobing error is 10µm for a 10mm radius part
form_2_max_amp = .0005     # worst case scenario for waviness error is <5µm for a 10mm radius part

init_measurement_pitch = 1 * (2 * np.pi / 360)


def fit_circle(x_list, y_list):
    # x_m = np.mean(x_list)
    # y_m = np.mean(y_list)

    def calc_r(xc, yc):
        """ calculate the distance of each 2D points from the center (xc, yc) """
        return np.sqrt((x_list - xc) ** 2 + (y_list - yc) ** 2)

    def f_2(c):
        """ calculate the algebraic distance between the data points and the mean circle centered at c=(xc, yc) """
        Ri = calc_r(*c)
        return Ri - Ri.mean()

    # center_estimate = x_m, y_m
    center_estimate = 0, 0
    # print(optimize.least_squares(f_2, center_estimate))
    center_2 = optimize.least_squares(f_2, center_estimate)['x']

    xc_2, yc_2 = center_2
    Ri_2 = calc_r(*center_2)
    R_2 = Ri_2.mean()
    residu_2 = sum((Ri_2 - R_2) ** 2)

    fit_results = {'x': xc_2, 'y': yc_2, 'r': R_2, 'residual': residu_2}
    return fit_results
    # return xc_2, yc_2, R_2


def calc_stuff(nom_diameter=init_nom_radius,
               form_1_amplitude=init_form_1_amp,
               form_2_amplitude=init_form_2_amp,
               form_1_wavelength=init_form_1_wave,
               form_2_wavelength=init_form_2_wave,
               start_angle=init_start_angle,
               measured_degrees=init_measured_degrees,
               measurement_pitch=init_measurement_pitch):
    #
    calc_results = {'circle_actual_xs': [],
                    'circle_actual_ys': [],
                    'arc_actual_xs': [],
                    'arc_actual_ys': [],
                    'arc_nominal_xs': [],
                    'arc_nominal_ys': [],
                    'arc_actual_residual': 1}

    for theta in np.arange(0, 2 * np.pi, 2 * np.pi / 360):
        r = nom_diameter + form_1_amplitude * np.sin(form_1_wavelength * theta)
        r = r * (nom_diameter + form_2_amplitude * np.sin(form_2_wavelength * theta))
        x = r * np.cos(theta)
        y = r * np.sin(theta)
        calc_results['circle_actual_xs'].append(x)
        calc_results['circle_actual_ys'].append(y)

    for meas_theta in np.arange(start_angle, start_angle + measured_degrees, measurement_pitch):
        r = nom_diameter + form_1_amplitude * np.sin(form_1_wavelength * meas_theta)
        r = r * (nom_diameter + form_2_amplitude * np.sin(form_2_wavelength * meas_theta))

        x = r * np.cos(meas_theta)
        y = r * np.sin(meas_theta)
        calc_results['arc_actual_xs'].append(x)
        calc_results['arc_actual_ys'].append(y)
        calc_results['arc_nominal_xs'].append(init_nom_radius * np.cos(meas_theta))
        calc_results['arc_nominal_ys'].append(init_nom_radius * np.sin(meas_theta))

    circle_actual_results = fit_circle(calc_results['circle_actual_xs'], calc_results['circle_actual_ys'])
    circle_actual_x = circle_actual_results['x']
    circle_actual_y = circle_actual_results['y']
    circle_actual_r = circle_actual_results['r']

    arc_actual_results = fit_circle(calc_results['arc_actual_xs'], calc_results['arc_actual_ys'])
    arc_actual_x = arc_actual_results['x']
    arc_actual_y = arc_actual_results['y']
    arc_actual_r = arc_actual_results['r']

    calc_results['arc_actual_residual'] = arc_actual_results['residual']

    nominal_full_circle = plt.Circle((circle_actual_x, circle_actual_y),
                                     circle_actual_r,
                                     color='g',
                                     fill=False,
                                     label='nominal_circle')

    calc_results['nominal_full_circle'] = nominal_full_circle

    arc_calc_circle = plt.Circle((arc_actual_x, arc_actual_y),
                                 arc_actual_r,
                                 color='y',
                                 fill=False,
                                 label='calculated_circle')

    calc_results['arc_calc_circle'] = arc_calc_circle

    return calc_results


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
                                  valinit=init_form_1_wave,
                                  valstep=1,
                                  orientation='vertical')

form_1_amplitude_slider = Slider(label='form_1\namplitude',
                                 ax=ax_form_1_amp,
                                 valmin=0,
                                 valmax=0.1,
                                 valinit=init_form_1_amp,
                                 # valstep=1,
                                 orientation='vertical')

form_2_wavelength_slider = Slider(label='form_2\nwavelengths',
                                  ax=ax_form_2_wavelength,
                                  valmin=1,
                                  valmax=180,
                                  valinit=init_form_2_wave,
                                  valstep=1,
                                  orientation='vertical')

form_2_amplitude_slider = Slider(label='form_2\namplitude',
                                 ax=ax_form_2_amp,
                                 valmin=0,
                                 valmax=0.1,
                                 valinit=init_form_2_amp,
                                 # valstep=1,
                                 orientation='vertical')


init_calc_results = calc_stuff()


ax.plot(init_calc_results['circle_actual_xs'], init_calc_results['circle_actual_ys'], label='actual_circle', color='r')
ax.plot(init_calc_results['arc_nominal_xs'], init_calc_results['arc_nominal_ys'], label='nominal_radius', color='b')
ax.plot(init_calc_results['arc_actual_xs'], init_calc_results['arc_actual_ys'], label='measured_points', color='k')
ax.add_patch(init_calc_results['nominal_full_circle'])
ax.add_patch(init_calc_results['arc_calc_circle'])

# ax.plot(full_xs, full_ys, label='actual_circle', color='r')
# ax.plot(nominal_xs, nominal_ys, label='nominal_radius', color='b')
# ax.plot(actual_xs, actual_ys, label='measured_points', color='k')
# ax.add_patch(nominal_full_circle)
# ax.add_patch(calc_circle)
# ax.legend(loc='upper left')


def update(val):
    ax.cla()
    start_angle = start_angle_slider.val * (2 * np.pi / 360)
    degrees_arc = dgrs_arc_slider.val * (2 * np.pi / 360)
    form_1_wavelengths = int(form_1_wavelength_slider.val)
    form_1_amp = form_1_amplitude_slider.val
    form_2_wavelengths = int(form_2_wavelength_slider.val)
    form_2_amp = form_2_amplitude_slider.val

    # full_xs, full_ys, actual_xs, actual_ys, nominal_xs, nominal_ys, nominal_full_circle, calc_circle = calc_stuff(
    #     start_angle=start_angle,
    #     measured_degrees=degrees_arc,
    #     form_1_amplitude=form_1_amp,
    #     form_1_wavelength=form_1_wavelengths,
    #     form_2_amplitude=form_2_amp,
    #     form_2_wavelength=form_2_wavelengths)

    new_calc_results = calc_stuff(start_angle=start_angle,
                                  measured_degrees=degrees_arc,
                                  form_1_amplitude=form_1_amp,
                                  form_1_wavelength=form_1_wavelengths,
                                  form_2_amplitude=form_2_amp,
                                  form_2_wavelength=form_2_wavelengths,
                                  measurement_pitch=init_measurement_pitch)

    ax.plot(new_calc_results['circle_actual_xs'], new_calc_results['circle_actual_ys'],
            label='actual_circle', color='r')
    ax.plot(new_calc_results['arc_nominal_xs'], new_calc_results['arc_nominal_ys'],
            label='nominal_radius', color='b')
    ax.plot(new_calc_results['arc_actual_xs'], new_calc_results['arc_actual_ys'],
            label='measured_points', color='k')
    ax.add_patch(new_calc_results['nominal_full_circle'])
    ax.add_patch(new_calc_results['arc_calc_circle'])

    # ax.plot(full_xs, full_ys, label='actual_circle', color='r')
    # ax.plot(nominal_xs, nominal_ys, label='nominal_radius', color='b')
    # ax.plot(actual_xs, actual_ys, label='measured_points', color='k')
    # ax.add_patch(nominal_full_circle)
    # ax.add_patch(calc_circle)
    ax.legend(loc='upper left')
    fig.canvas.draw()


start_angle_slider.on_changed(update)
dgrs_arc_slider.on_changed(update)
form_1_wavelength_slider.on_changed(update)
form_1_amplitude_slider.on_changed(update)
form_2_wavelength_slider.on_changed(update)
form_2_amplitude_slider.on_changed(update)

plt.show()
