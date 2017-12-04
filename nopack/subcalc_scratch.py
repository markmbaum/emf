#import cProfile
#import pstats

import emf.subcalc as sc

mod = sc.load_towers(r"C:\Users\mab0328\Documents\code\emf\working-files\subcalc\TEST1.INP", True)
mod.N = 25e2
mod.calculate()

#cProfile.run(code, filename='profile')
#pstats.Stats('profile').strip_dirs().sort_stats('time').print_stats(50)
