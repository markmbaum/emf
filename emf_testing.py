import emf
import cProfile
import pstats

#emf.run('practice_xcs.xlsx', path = 'run-dest/')

sb = emf.load_template('practice_xcs.xlsx')

xc = sb['32P']

cProfile.run('emf.optimize_phasing(xc)',
filename = 'profile')

pstats.Stats('profile').strip_dirs().sort_stats('time').print_stats(50)
