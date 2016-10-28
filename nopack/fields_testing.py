import cProfile
import pstats

import emf.fields as fld

sb = fld.load_template(r"P:\MBaum\Programming\Python\python_code\emf\working_files\practice_xcs.xlsx")

fld.plot_xs(sb['raise1'])
fld.show()

#cProfile.run("""
#xc = fields.load_template('../working_files/practice_xcs.xlsx').sample()
#res, opt = fields.optimize_phasing(xc, 'all')
#""",
#filename = 'profile')

#pstats.Stats('profile').strip_dirs().sort_stats('time').print_stats(50)
