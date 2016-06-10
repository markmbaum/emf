import emf
import cProfile
import pstats

#sb = emf.run('practice_xcs.xlsx', path = 'run-dest/')
#sb = emf.load_template('practice_xcs.xlsx')
#opt = emf.optimize_phasing(sb['HL_P'])
#opt.to_excel('optimize_phasing_test.xlsx', index = False)

cProfile.run("""sb = emf.load_template('practice_xcs.xlsx')""",
    filename = 'profile')
pstats.Stats('profile').strip_dirs().sort_stats('time').print_stats(50)
