import cProfile
import pstats
import sys

from emf import subcalc as sc

mod = sc.load_model(r"P:\MBaum\Programming\Python\python_code\emf\working_files\REF_GRID1.REF",
		r"P:\MBaum\Programming\Python\python_code\emf\working_files\footprints1.csv")

x,y,b = mod.path((10,20,30,400), (10,20,30,400), n=100)

print b
