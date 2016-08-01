import emf
import cProfile
import pstats


cProfile.run("""
xc = emf.load_template('practice_xcs').sample()
res, opt = emf.optimize_phasing(xc, 'all', path = 'dump')
""",
filename = 'profile')

pstats.Stats('profile').strip_dirs().sort_stats('time').print_stats(50)
