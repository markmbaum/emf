import cProfile
import pstats
import sys

sys.path.append("/Volumes/BAUMPRIMARY/Code/Python/")
sys.path.append("/Volumes/BAUMPRIMARY/Code/Python/emf/")

import emf

cProfile.run("""
xc = emf.load_template('/Volumes/BAUMPRIMARY/Code/Python/emf/working_files/practice_xcs.xlsx').sample()
res, opt = emf.optimize_phasing(xc, 'all')
""",
filename = 'profile')

pstats.Stats('profile').strip_dirs().sort_stats('time').print_stats(50)
