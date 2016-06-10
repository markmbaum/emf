import emf
import cProfile
import pstats

#emf.run('practice_xcs.xlsx', path = 'run-dest/')

sb = emf.load_template('practice_xcs.xlsx')

xc = sb['32P']
print(xc.y)
print(emf.target_fields(xc, None, 10., None, .05, range((len(xc.hot))), range(len(xc.gnd))))

"""cProfile.run(,
filename = 'profile')

pstats.Stats('profile').strip_dirs().sort_stats('time').print_stats(50)"""
