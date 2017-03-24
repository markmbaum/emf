"""
The emf package is a container for two subpackages. The subpackages, which are called fields and subcalc, are for modeling electric and magnetic fields near power lines. The fields subpackage is for 2D models. It predicts electric and magnetic fields near sets of parallel power lines. The subcalc package is for 3D models, predicting the magnetic fields generated by complex, nonparallel collections of current carrying wire segments. Both subpackages are freestanding, meaning they can be used to build, compute, and analyze models all by themselves without any other programs. They rely heavily on numpy, pandas, and matplotlib, but the EMF calculations do not rely on algorithms from any external libraries or legacy code.

Both fields and subcalc were originally meant to supplement older modeling programs and eventually became complete replacements of them. Thus, both packages are styled after the older programs and still contain a fair amount of functionality for working with those programs. Both packages also have specifically formatted excel/csv files for storing input data and results. Although these template files are never absolutely necessary, they make the packages much more usable for those who only know a little bit of Python and they make larger-scale modeling projects more repeatable and better documented. Neither package has a fancy UI, just scripts and (optionally) template files.
"""

import os
import copy
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
from mpl_toolkits.mplot3d import Axes3D as _Axes3D
from matplotlib.patches import Rectangle as _Rectangle
from pkg_resources import resource_filename as _resource_filename

import fields
import subcalc

__all__ = ['fields', 'subcalc']
