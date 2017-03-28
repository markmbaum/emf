#import cProfile
#import pstats

import emf.subcalc as sc

res = sc.load_results(r"G:\Projects\215106_NeedhamEMF\Models\second-round-ir\subcalc\R18.REF")

sc.plot_contour(res)
sc.show()

#cProfile.run(code, filename='profile')
#pstats.Stats('profile').strip_dirs().sort_stats('time').print_stats(50)
