import emf
import cProfile
import pstats

#emf.run('practice_xcs.xlsx', path = 'run-dest/')
sb = emf.load_template('practice_xcs')

fn = r"P:\MBaum\Programming\Python\python_code\FLD\XC-comparisons\32P.DAT"

sb['32P'].compare_DAT(fn, save = True)

"""cProfile.run('emf.optimize_phasing(xc)',
filename = 'profile')

pstats.Stats('profile').strip_dirs().sort_stats('time').print_stats(50)
"""
