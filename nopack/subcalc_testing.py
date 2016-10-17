import cProfile
import pstats
import sys

from emf import subcalc as sc

mod = sc.load_model(r"P:\MBaum\Programming\Python\python_code\emf\working_files\REF_GRID1.REF",
		r"P:\MBaum\Programming\Python\python_code\emf\working_files\footprints1.csv")

fig, ax, CS = sc.plot_contour(mod, scale='log', levels=[.1,.5,1,5,10,25,50])
sc.show()
sc.close()

mod = sc.load_model(r"P:\MBaum\Programming\Python\python_code\emf\working_files\REF_GRID2.xlsx")
mod.Bkey = 'Bres'
fig, ax, QM, cbar = sc.plot_pcolormesh(mod)
sc.show()
