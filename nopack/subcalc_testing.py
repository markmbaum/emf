#import cProfile
#import pstats

import emf.subcalc as sc

mod = sc.read_INP(r"P:\MBaum\Programming\Python\python_code\emf\working_files\subcalc\SUB211P.INP", return_model=True)

mod.xlim = 350, 1200
mod.ylim = 300, 850
mod.spacing = 1.25

sc.plot_cross_sections(
    mod.calculate(),
    [[(400, 700), (1000, 800)], [(400, 400), (700, 350), (1150, 550)]],
    n=250,
    max_fig_height=10,
    max_fig_width=16,
    x_labeling='location',
    levels=[0.1, 0.5, 1, 5, 10, 25, 50, 100],
    label_max=False,
    save=True)

sc.show()

#cProfile.run(code, filename='profile')
#pstats.Stats('profile').strip_dirs().sort_stats('time').print_stats(50)
