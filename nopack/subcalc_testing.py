import cProfile
import pstats
import sys

sys.path.append("/Volumes/BAUMPRIMARY/Code/Python/")
sys.path.append("/Volumes/BAUMPRIMARY/Code/Python/emf/")

from emf import subcalc

fn = '../working_files/REF_GRID.xlsx'

mod = subcalc.load_model(fn)

mod.north_angle = 5

subcalc.plot_contours(mod, levels=[0.1,0.5,1,5,10,25,50])

subcalc.show()
