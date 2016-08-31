import cProfile
import pstats
import sys

sys.path.append("/Volumes/BAUMPRIMARY/Code/Python/")
sys.path.append("/Volumes/BAUMPRIMARY/Code/Python/emf/")

from emf import fields

cProfile.run("""
xc = fields.load_template('../working_files/practice_xcs.xlsx').sample()
res, opt = fields.optimize_phasing(xc, 'all')
""",
filename = 'profile')

pstats.Stats('profile').strip_dirs().sort_stats('time').print_stats(50)
