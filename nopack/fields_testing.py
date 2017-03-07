import cProfile
import pstats

from emf import fields as fld

print fld.read_FLDs(r'G:\Projects\215106_NeedhamEMF\Models\initial-modeling\FIELDS input and output files')

#cProfile.run("""
#xc = fields.load_template('../working_files/practice_xcs.xlsx').sample()
#res, opt = fields.optimize_phasing(xc, 'all')
#""",
#filename = 'profile')

#pstats.Stats('profile').strip_dirs().sort_stats('time').print_stats(50)
