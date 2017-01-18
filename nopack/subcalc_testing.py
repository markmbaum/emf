#import cProfile
#import pstats

import numpy as np

import emf.subcalc as sc

res = sc.load_results(r"P:\MBaum\Programming\Python\python_code\emf\working_files\subcalc\REF_GRID1.xlsx")

#cProfile.run(code, filename='profile')
#pstats.Stats('profile').strip_dirs().sort_stats('time').print_stats(50)
