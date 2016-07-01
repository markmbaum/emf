import emf
import cProfile
import pstats

sb = emf.load_template('practice_xcs.xlsx')

xc = sb['HL_P']

emf.optimize_phasing(xc, 'all', path = 'ophasing/')

"""
cProfile.run('emf.optimize_phasing(xc)',
filename = 'profile')

pstats.Stats('profile').strip_dirs().sort_stats('time').print_stats(50)
"""
