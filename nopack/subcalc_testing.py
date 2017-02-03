#import cProfile
#import pstats

import pandas as pd

import emf.subcalc as sc

mod = sc.load_towers(r"G:\Projects\215106_NeedhamEMF\Models\first-round-ir\subcalc\baker-st-substation\baker-st-towers.xlsx", True, sheet='enorm')

#get footprint df
fps = pd.read_csv(r"G:\Projects\215106_NeedhamEMF\Models\first-round-ir\subcalc\baker-st-substation\baker-st-fence.csv")

#select the fenceline points for plots
pts = zip(fps['X'].values[5:][::-1], fps['Y'].values[5:][::-1])

sc.plot_path(mod, pts, cmap='magma')

sc.show()

#cProfile.run(code, filename='profile')
#pstats.Stats('profile').strip_dirs().sort_stats('time').print_stats(50)
