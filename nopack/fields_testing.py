import cProfile
import pstats

import emf.fields as fld

c = fld.Conductor('1a', dict(x=10, y=23, V=345, I=200))
print c
print fld.Conductor('1b', dict(x=40), c)

#cProfile.run("""
#xc = fields.load_template('../working_files/practice_xcs.xlsx').sample()
#res, opt = fields.optimize_phasing(xc, 'all')
#""",
#filename = 'profile')

#pstats.Stats('profile').strip_dirs().sort_stats('time').print_stats(50)
