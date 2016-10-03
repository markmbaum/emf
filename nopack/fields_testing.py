import cProfile
import pstats

from emf import fields as fld

sb = fld.load_template('../working_files/practice_xcs')

fld.plot_ROW_values(sb, path='test')

#cProfile.run("""
#xc = fields.load_template('../working_files/practice_xcs.xlsx').sample()
#res, opt = fields.optimize_phasing(xc, 'all')
#""",
#filename = 'profile')

#pstats.Stats('profile').strip_dirs().sort_stats('time').print_stats(50)
