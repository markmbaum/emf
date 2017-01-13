#import cProfile
#import pstats

import numpy as np

import emf.subcalc as sc

A = [(20, 30, 10),
    (20, 20, 20),
    (20, 10, 10),
    (50, 60, 10),
    (50, 50, 20),
    (50, 40, 10)]

B = [(50, 60, 10),
    (50, 50, 20),
    (50, 40, 10),
    (90, 60, 30),
    (90, 50, 40),
    (90, 40, 50)]

Ph = [0., -120., 120., 0., -120., 120.]

A = [np.array(a, dtype=float) for a in A]
B = [np.array(b, dtype=float) for b in B]
I = 100.

x = np.linspace(0, 100, 500)
y = np.linspace(0, 100, 500)
z = 3.28

Ph_x, Ph_y, Ph_z = sc.B_field_segment(A[0], B[0], I, Ph[0], x, y, z)
for (a,b,ph) in zip(A[1:], B[1:], Ph[1:]):
    Ph = sc.B_field_segment(a, b, I, ph, x, y, z)
    Ph_x += Ph[0]
    Ph_y += Ph[1]
    Ph_z += Ph[2]

Ph_x, Ph_y, Ph_z, X, Y = sc.grid_segment_results(Ph_x, Ph_y, Ph_z, x, y)

Bx, By, Bz, Bres, Bmax = sc.phasors_to_magnitudes(Ph_x, Ph_y, Ph_z)

mod = sc.Results(dict(X=X, Y=Y, Bx=Bx, By=By, Bz=Bz, Bres=Bres, Bmax=Bmax))

mod.Bkey = 'Bx'

sc.plot_contour(mod, max_fig_height=5, save=True)

#cProfile.run(code, filename='profile')
#pstats.Stats('profile').strip_dirs().sort_stats('time').print_stats(50)
