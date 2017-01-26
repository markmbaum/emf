"""
The emf package is a container for two subpackages. The subpackages, which are called fields and subcalc, are for modeling electric and magnetic fields near power lines. The fields subpackage is for 2D models and the subcalc package is for 3D models. Both subpackages can be used to build, compute, and analyze models. They rely heavily on numpy, pandas, and matplotlib, but the EMF calculations do not rely on any external libraries or legacy code. The 2D modeling calculations are done in Python (with numpy) and the 3D calculations are done with a combination of Python and C. Brief explanations of the packages are below and links to each subpackage's documentation are on the left.

1. emf.fields was originally a small stand-alone package that streamlined the use of an old electromagnetic field (EMF) modeling program called FIELDS. That old program predicts electric and magnetic fields near parallel sets of power lines by assuming the conductors are infinitely long and computing the fields along a transect perpendicular to the power lines (a cross section model). The old FIELDS program is very difficult to use, so it made sense to transplant as much of the modeling process as possible into other programs. So, emf.fields was initially an extension of the old FIELDS program that made it easier to work with the modeling files and do things like plot the results. In its first versions, emf.fields was used to read power line information from specially formatted excel templates, create input files that could be used by FIELDS to compute EMF results, then manage the FIELDS results. The package did everything except perform the actual EMF calculations. In its current versions, emf.fields still contains functions streamlining the use of FIELDS and it performs all of the necessary calculations on its own. It also has a significant amount of fuctionality that FIELDS lacks, like functions to optimize the phasing arrangement of a cross section model, compute conductor heights needed to lower fields to target levels, and create a handful of different plots. For someone that knows a little Python, emf.fields is signifantly more flexible and more powerful than the old FIELDS program.

2. Like emf.fields, emf.subcalc was originally intended to supplement/extend the use of another modeling program. In this case, the other modeling program is called SUBCALC (developed by Enertech, sponsored by EPRI. SUBCALC predicts magnetic fields over a fixed-height 2 dimensional grid and can model many non-parallel segments of power lines. It is generally used to model magnetic fields near electrical substations, where circuit arrangements can be complex. Just like emf.fields, emf.subcalc is now a complete replacement of the program it was originally intended to supplement. Although it provides no fancy GUI, emf.subcalc is capable of calculating the magnetic fields produced by complex, three dimensional arrangements of current carrying wire segments.
"""

import os
import copy
import glob
import ctypes
import shutil
import datetime
import textwrap
import itertools
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.interpolate import interpn as _interpn

import fields
import subcalc

__all__ = ['fields', 'subcalc']
