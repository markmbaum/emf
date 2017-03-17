#import cProfile
#import pstats

import pandas as pd

import emf.subcalc as sc

res = sc.load_towers('../working_files/subcalc/TEST1.INP', True).calculate()

print res.info

#cProfile.run(code, filename='profile')
#pstats.Stats('profile').strip_dirs().sort_stats('time').print_stats(50)
