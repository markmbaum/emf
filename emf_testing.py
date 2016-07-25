import emf
import cProfile
import pstats

emf.convert_DAT_crawl('*', bundle = True)

"""
cProfile.run('emf.optimize_phasing(xc)',
filename = 'profile')

pstats.Stats('profile').strip_dirs().sort_stats('time').print_stats(50)
"""
