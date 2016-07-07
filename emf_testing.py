import emf
import cProfile
import pstats

sb = emf.load_template(r"G:\Projects\215030_Millbury EMF NGrid\Models\FIELDS\Updated FIELDS 070716\Millbury_Average_Loading.xlsx")

xc = sb.i(0)

emf.target_fields(xc,.01,.01,.01,.01,'all','all',path = r"G:\Projects\215030_Millbury EMF NGrid\Models\FIELDS\Updated FIELDS 070716\New folder")
"""
cProfile.run('emf.optimize_phasing(xc)',
filename = 'profile')

pstats.Stats('profile').strip_dirs().sort_stats('time').print_stats(50)
"""
