#import cProfile
#import pstats

from emf import subcalc as sc

mod = sc.load_model(r"G:\Projects\215106_NeedhamEMF\Models\first-round-ir\subcalc\baker-st-substation\PNORM.REF", Bkey='Bres')

points_a = ((50,100), (200,30), (300,40), (400,200))
points_b = ((102,102), (22,302), (302,402), (400,20))


sc.plot_cross_sections(mod, (points_a, points_b), max_fig_width=16, x_labeling='location', xs_label_size=16, n=500)

sc.show()
