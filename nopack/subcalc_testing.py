#import cProfile
#import pstats

import emf.subcalc as sc

res = sc.load_results(r"P:\MBaum\Programming\Python\python_code\emf\working_files\subcalc\TEST2.REF")

N = res.N

resres = res.resample(N=int(N/250.0))

print res.N, resres.N

sc.plot_contour(res)
sc.plot_contour(resres)

sc.show()

#cProfile.run(code, filename='profile')
#pstats.Stats('profile').strip_dirs().sort_stats('time').print_stats(50)
