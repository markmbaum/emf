import cProfile
import pstats
import sys

sys.path.append("/Volumes/BAUMPRIMARY/Code/Python/")
sys.path.append("/Volumes/BAUMPRIMARY/Code/Python/emf/")

from emf import subcalc

fn_data = '../working_files/REF_GRID'
fn_footprints = '../working_files/footprints'

mod = subcalc.load_model(fn_data, fn_footprints)

fig,ax,CS = subcalc.contour_plot(mod)

ax.grid(b = True)

subcalc.show()
