import emf
import cProfile
import pstats

sb = emf.load_template(r"\\camfs\G_Drive\Projects\215030_Millbury EMF NGrid\Models\FIELDS\Updated FIELDS 070716\Millbury_Average_Loading.xlsx")

xc = sb.i(-1)

res,opt = emf.optimize_phasing(xc, 'all', save = True)

"""
cProfile.run('emf.optimize_phasing(xc)',
filename = 'profile')

pstats.Stats('profile').strip_dirs().sort_stats('time').print_stats(50)
"""
