import numpy as np
from scipy import optimize
import csv
import progressbar

widgets = [' [', progressbar.Timer(format='elapsed time: %(elapsed)s'), '] ', progressbar.Bar('*'),
           ' (', progressbar.AdaptiveETA(), ') ', ]

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

    # calc_results = {'circle_actual_xs': [],
    #                 'circle_actual_ys': [],
    #                 'arc_actual_xs': [],
    #                 'arc_actual_ys': [],
    #                 'arc_nominal_xs': [],
    #                 'arc_nominal_ys': [],
    #                 'arc_actual_residual': 1}

    circle_actual_xs = []
    circle_actual_ys = []
    arc_actual_xs = []
    arc_actual_ys = []

    calc_results = {'circle_actual_x': 0,
                    'circle_actual_y': 0,
                    'circle_actual_r': 1,
                    'arc_actual_x': 0,
                    'arc_actual_y': 0,
                    'arc_actual_r': 0,
                    'arc_actual_residual': 1,
                    'position_error': 0,
                    'size_error': 0}

    # for theta in np.arange(0, 2 * np.pi, 2 * np.pi / 360):
    #     r = nom_diameter + form_1_amplitude * np.sin(form_1_wavelength * theta)
    #     r = r * (nom_diameter + form_2_amplitude * np.sin(form_2_wavelength * theta))
    #     x = r * np.cos(theta)
    #     y = r * np.sin(theta)
    #     circle_actual_xs.append(x)
    #     circle_actual_ys.append(y)

    for meas_theta in np.arange(start_angle, start_angle + measured_degrees, measurement_pitch):
        r = nom_diameter + form_1_amplitude * np.sin(form_1_wavelength * meas_theta)
        r = r * (nom_diameter + form_2_amplitude * np.sin(form_2_wavelength * meas_theta))

        x = r * np.cos(meas_theta)
        y = r * np.sin(meas_theta)
        arc_actual_xs.append(x)
        arc_actual_ys.append(y)

    # circle_actual_results = fit_circle(circle_actual_xs, circle_actual_ys)
    # calc_results['circle_actual_x'] = circle_actual_results['x']
    # calc_results['circle_actual_y'] = circle_actual_results['y']
    # calc_results['circle_actual_r'] = circle_actual_results['r']

    arc_actual_results = fit_circle(arc_actual_xs, arc_actual_ys)
    calc_results['arc_actual_x'] = arc_actual_results['x']
    calc_results['arc_actual_y'] = arc_actual_results['y']
    calc_results['arc_actual_r'] = arc_actual_results['r']

    calc_results['arc_actual_residual'] = arc_actual_results['residual']

    position_error = np.sqrt((calc_results['circle_actual_x'] - calc_results['arc_actual_x']) ** 2 +
                             (calc_results['circle_actual_y'] - calc_results['arc_actual_y']) ** 2)
    calc_results['position_error'] = position_error

    size_error = calc_results['circle_actual_r'] - calc_results['arc_actual_r']
    calc_results['size_error'] = size_error

    return calc_results

# results = []

output_file = r'C:\Users\benhartd\PycharmProjects\random\radius_analysis_output.csv'

header = ['start_angle',
          'degrees_of_arc',
          'pitch',
          'form_1_amp',
          'form_1_wave',
          'circle_actual_x',
          'circle_actual_y',
          'circle_actual_r',
          'arc_actual_x',
          'arc_actual_y',
          'arc_actual_r',
          'arc_actual_residual',
          'position_error',
          'size_error']


bar_max = 18227430
print(bar_max)

bar = progressbar.ProgressBar(max_value=bar_max, widgets=widgets).start()

i = 0

with open(output_file, 'w') as csvfile:
    csvwriter = csv.writer(csvfile)

    for form_1_wave in range(2, 5):
        max_start_angle = 360 / form_1_wave * (2 * np.pi / 360)
        for start_angle in np.arange(0, max_start_angle, 2 * np.pi / 360):
            for degrees_of_arc in np.arange(10 * (2 * np.pi / 360), 179 * (2 * np.pi / 360), 2 * np.pi / 360):
                max_pitch = degrees_of_arc / 3
                for pitch in np.arange(1 * (2 * np.pi / 360), max_pitch, 1 * (2 * np.pi / 360)):
                    # print(start_angle, degrees_of_arc, pitch)
                    for form_1_amp in np.arange(0.0001, 0.001, 0.0001):
                        i += 1
                        bar.update(i)
                        result = [start_angle, degrees_of_arc, pitch, form_1_amp, form_1_wave]
                        calc_results = calc_stuff(start_angle=start_angle,
                                                  measured_degrees=degrees_of_arc,
                                                  measurement_pitch=pitch,
                                                  form_1_amplitude=form_1_amp,
                                                  form_1_wavelength=form_1_wave,
                                                  form_2_amplitude=0,
                                                  form_2_wavelength=1)

                        result.append(calc_results['circle_actual_x'])
                        result.append(calc_results['circle_actual_y'])
                        result.append(calc_results['circle_actual_r'])
                        result.append(calc_results['arc_actual_x'])
                        result.append(calc_results['arc_actual_y'])
                        result.append(calc_results['arc_actual_r'])
                        result.append(calc_results['arc_actual_residual'])
                        result.append(calc_results['position_error'])
                        result.append(calc_results['size_error'])

                        csvwriter.writerow(result)
print(i)
