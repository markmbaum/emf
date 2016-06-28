import emf
import cProfile
import pstats

emf.run(r"FIELDS_StonehamIRs.xlsx", path = r'FIELDS_StonhamIRs/')

"""cProfile.run('emf.optimize_phasing(xc)',
filename = 'profile')

pstats.Stats('profile').strip_dirs().sort_stats('time').print_stats(50)
"""
