import cProfile
import pstats
import sys

from emf import subcalc as sc

fn_grid = '../working_files/REF_GRID.REF'
fn_foot = '../working_files/footprints.csv'

sc.drop_template()
