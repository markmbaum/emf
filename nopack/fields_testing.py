import cProfile
import pstats

from emf import fields as fld

fn = r"G:\Projects\EMF Resources\fields-practice\MMB\practice-cross-section.xlsx"

fld.to_FLDs(fn)

#cProfile.run("""
#xc = fields.load_template('../working_files/practice_xcs.xlsx').sample()
#res, opt = fields.optimize_phasing(xc, 'all')
#""",
#filename = 'profile')

#pstats.Stats('profile').strip_dirs().sort_stats('time').print_stats(50)
