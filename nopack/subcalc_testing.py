import cProfile
import pstats
import sys

from emf import subcalc

fn_grid = '../working_files/REF_GRID.REF'
fn_foot = '../working_files/footprints.csv'

mod = subcalc.load_model(fn_grid, fn_foot)

mod.north_angle = 5

subcalc.plot_pcolormesh(mod, path='../docs/images/pcolormesh_plot')

subcalc.show()
