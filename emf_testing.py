import emf
import cProfile
import pstats

sb = emf.load_template('practice_xcs.xlsx')

cProfile.run("""
opt = emf.optimize_phasing(sb['HL_P'])
""",
filename = 'profile')

opt.to_excel('optimize_phasing_test2.xlsx')

pstats.Stats('profile').strip_dirs().sort_stats('time').print_stats(50)
