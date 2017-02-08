import cProfile
import pstats

from emf import fields as fld

c = fld.Conductor('z', [1,1,1,1,1,1,1,1])
xs = fld.CrossSection('x', [c])
xs.max_dist = 100
xs.lROW = -50
xs.rROW = 50
xs.step = 0.1
sb = fld.SectionBook('test', [xs])

fld.plot_max_fields(sb.i[0])
fld.show()

#cProfile.run("""
#xc = fields.load_template('../working_files/practice_xcs.xlsx').sample()
#res, opt = fields.optimize_phasing(xc, 'all')
#""",
#filename = 'profile')

#pstats.Stats('profile').strip_dirs().sort_stats('time').print_stats(50)
