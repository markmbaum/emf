#import cProfile
#import pstats

import emf.subcalc as sc

mod = sc.load_towers(r"G:\Projects\EMF_Resources\mmb-python-package\working_files\subcalc\TEST1.INP", True)

print mod.sample(1,1,range(1000))

#cProfile.run(code, filename='profile')
#pstats.Stats('profile').strip_dirs().sort_stats('time').print_stats(50)
