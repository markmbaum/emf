#import cProfile
#import pstats

import emf.subcalc as sc

mod = sc.Model(name='test',
    towers=[
        sc.Tower('line1', 0, 10, 10, -90, [-1,0,1], [10,20,10], [10,11,12], [0,120,240]),
        sc.Tower('line1', 1, 140, 200, 0, [-10,0,10], [20,10,20], [10,11,12], [0,120,240]),
        sc.Tower('line1', 2, 300, 250, -45, [-1,0,1], [20,20,20], [10,11,12], [0,120,240]),
        sc.Tower('line2', 0, 90, 300, 0, [-3,-5,-3], [10,20,30], [10,11,12], [0,120,240]),
        sc.Tower('line2', 1, 110, 250, 0, [-3,-5,-3], [10,20,30], [10,11,12], [0,120,240]),
    ],
    spacing=1
)

res = mod.calculate()

sc.plot_pcolormesh(res, save=True)

#cProfile.run(code, filename='profile')
#pstats.Stats('profile').strip_dirs().sort_stats('time').print_stats(50)
