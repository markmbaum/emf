"""The emf package is a container for two subpackages, emf.fields and emf.subcalc. It imports libraries used by the subpackages like numpy, pandas, and matplotlib. It also contains some custom functions and classes that are used by both emf.fields and emf.subcalc, but they are all private. See the documentation for emf.fields and emf.subcalc."""

import os
import copy
import glob
import ctypes
import shutil
import textwrap
import itertools
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.integrate import quad as _quad
from scipy.interpolate import interpn as _interpn

import fields
import subcalc

__all__ = ['fields', 'subcalc']
