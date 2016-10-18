import cProfile
import pstats

import emf.fields as fld

sb = fld.load_template(r"P:\MBaum\Programming\Python\python_code\emf\working_files\practice_xcs.xlsx")

xs = sb['32P']

h, adj = fld.target_fields(xs, ['1a', '1b', '1c', '2a', '2b', '2c'], 100, 0, 1, 0,
		rel_err=1.0e-3, save=True)

print h, adj.ROW_edge_max

figs = fld.plot_groups_at_ROW(adj, return_figs=True)


#cProfile.run("""
#xc = fields.load_template('../working_files/practice_xcs.xlsx').sample()
#res, opt = fields.optimize_phasing(xc, 'all')
#""",
#filename = 'profile')

#pstats.Stats('profile').strip_dirs().sort_stats('time').print_stats(50)
