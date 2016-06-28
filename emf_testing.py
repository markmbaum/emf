import emf
import cProfile
import pstats

#emf.run(r"G:\Projects\215008_EMF_Woburn_Wakefield\Models\FIELDS_Stoneham\FIELDS_StonehamIRs.xlsx")

emf.DAT_to_csv_crawl(r"G:\Projects\215008_EMF_Woburn_Wakefield\Models\FIELDS_Stoneham")

"""cProfile.run('emf.optimize_phasing(xc)',
filename = 'profile')

pstats.Stats('profile').strip_dirs().sort_stats('time').print_stats(50)
"""
