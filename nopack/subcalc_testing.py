import cProfile
import pstats

code = """
import numpy as np
from emf import subcalc as sc

#circuit 1 points
p1A_a, p1A_b = np.array([50, 60, 10.0]), np.array([100, 110, 20.0])
p1B_a, p1B_b = np.array([50, 50, 20.0]), np.array([100, 100, 30.0])
p1C_a, p1C_b = np.array([50, 40, 10.0]), np.array([100, 90, 20.0])

#circuit 2 points
p2A_a, p2A_b = np.array([190, 200, 10.0]), np.array([190, 50, 20.0])
p2B_a, p2B_b = np.array([200, 200, 20.0]), np.array([200, 50, 30.0])
p2C_a, p2C_b = np.array([210, 200, 10.0]), np.array([210, 50, 20.0])

#samples
x = np.arange(251, dtype=float)
y = np.arange(251, dtype=float)
z = np.array([3.28])

#circuit 1 phasors
B_1A = sc.B_field_segment(p1A_a, p1A_b, 100, 0, x, y, z)
B_1B = sc.B_field_segment(p1B_a, p1B_b, 100, -120, x, y, z)
B_1C = sc.B_field_segment(p1C_a, p1C_b, 100, 120, x, y, z)

#circuit 2 phasors
B_2A = sc.B_field_segment(p2A_a, p2A_b, 100, 120, x, y, z)
B_2B = sc.B_field_segment(p2B_a, p2B_b, 100, -120, x, y, z)
B_2C = sc.B_field_segment(p2C_a, p2C_b, 100, 0, x, y, z)

#combine phasors
Ph_x = B_1A[0] + B_1B[0] + B_1C[0] + B_2A[0] + B_2B[0] + B_2C[0]
Ph_y = B_1A[1] + B_1B[1] + B_1C[1] + B_2A[1] + B_2B[1] + B_2C[1]
Ph_z = B_1A[2] + B_1B[2] + B_1C[2] + B_2A[2] + B_2B[2] + B_2C[2]

#calc magnitudes
Bx, By, Bz, Bres, Bmax = sc.phasors_to_magnitudes(Ph_x, Ph_y, Ph_z)

X, Y = np.meshgrid(x, y, indexing='xy')
Y = Y[::-1,:]

mod = sc.Model(dict(X=X, Y=Y,
        Bx=Bx[:,:,0],
        By=By[:,:,0],
        Bz=Bz[:,:,0],
        Bres=Bres[:,:,0],
        Bmax=Bmax[:,:,0]))

print mod.B"""

cProfile.run(code, filename='profile')
pstats.Stats('profile').strip_dirs().sort_stats('time').print_stats(50)
