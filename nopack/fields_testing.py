import cProfile
import pstats

from emf import fields as fld

fn2 = r"G:\Projects\EMF_Resources\mmb-python-package\working_files\fields\enviro\enviro-files\Envsmpl3.o02"

fn1 = r"G:\Projects\EMF_Resources\mmb-python-package\working_files\fields\enviro\enviro-files\Envsmpl3.o01"

pan, figs = fld.enviro.compare_o01_o02(fn1, fn2)

fld.show()

#cProfile.run("""
#xc = fields.load_template('../working_files/practice_xcs.xlsx').sample()
#res, opt = fields.optimize_phasing(xc, 'all')
#""",
#filename = 'profile')

#pstats.Stats('profile').strip_dirs().sort_stats('time').print_stats(50)
