#import cProfile
#import pstats

import numpy as np
import matplotlib.pyplot as plt

import emf.subcalc as sc

Nx, Ny = 225, 100
x = np.linspace(0, 100, Nx)
y = np.linspace(0, 100, Ny)
a = np.array([20, 50., 10.], dtype=float)
b = np.array([90, 20, 10.], dtype=float)

Ph_x_1, Ph_y_1, Ph_z_1 = sc.B_field_segment(a, b, 100., 0., x, y, 0.)
Ph_x_2, Ph_y_2, Ph_z_2 = sc.B_field_segment(a+5, b+10, 100., -0., x, y, 0.)
Ph_x_3, Ph_y_3, Ph_z_3 = sc.B_field_segment(a+5, b+10, 100., 0., x, y, 0.)
Ph_x = Ph_x_1 + Ph_x_2 + Ph_x_3
Ph_y = Ph_y_1 + Ph_y_2 + Ph_y_3
Ph_z = Ph_z_1 + Ph_z_2 + Ph_z_3

Ph_x, Ph_y, Ph_z, X, Y = sc.grid_segment_results(Ph_x, Ph_y, Ph_z, x, y)

B = sc.phasors_to_magnitudes(Ph_x, Ph_y, Ph_z)

mod = sc.Model(X, Y, B[-1])

sc.plot_cross_sections(mod, ([(10,10), (10,45), (60,50)],), path='../')

#cProfile.run(code, filename='profile')
#pstats.Stats('profile').strip_dirs().sort_stats('time').print_stats(50)
