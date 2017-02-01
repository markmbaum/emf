#import cProfile
#import pstats

import numpy as np

import emf.subcalc as sc

mod = sc.load_towers(r"P:\MBaum\Programming\Python\python_code\emf\docs\notebooks\subcalc\towers.csv", True)

mod.spacing = 5

df = mod.sample(1, np.linspace(10,1000,25), np.linspace(1, 100, 25))
#df = mod.sample(1, 2, 3)

print df

#cProfile.run(code, filename='profile')
#pstats.Stats('profile').strip_dirs().sort_stats('time').print_stats(50)
