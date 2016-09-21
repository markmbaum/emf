import cProfile
import pstats
import sys

from emf import subcalc

fn_grid = '../working_files/REF_GRID.REF'
fn_foot = '../working_files/footprints.csv'

mod = subcalc.load_model(fn_grid, fn_foot)

mod.north_angle = 5

subcalc.plot_contour(mod, scale='log', levels=[.1,.5,1,5,10,25,50],
        path='../docs/images/contour_plot_log')

subcalc.show()
