# `emf`

![biot-savart](docs/img/both-equations.png)

The `emf` package is a container for two subpackages (documentation accessible [here](http://mbaum1122.github.io/emf/)):

1. `emf.fields` was originally a small stand-alone package that streamlined the use of an old electromagnetic field (EMF) modeling program called FIELDS. That old program predicts electric and magnetic fields near parallel sets of power lines by assuming the conductors are infinitely long and computing the fields along a transect perpendicular to the power lines (a cross section model). The old FIELDS program is very difficult to use (details on why below), so it made sense to transplant as much of the modeling process as possible into other programs. So, `emf.fields` was initially an extension of the old FIELDS program that made it easier to work with the modeling files and do things like plot the results. In its first versions, `emf.fields` was used to read power line information from specially formatted excel templates, create input files that could be used by FIELDS to compute EMF results, then manage the FIELDS results. The package did everything except perform the actual EMF calculations. In its current versions, `emf.fields` still contains functions streamlining the use of FIELDS *and* it performs all of the necessary calculations on its own. It also has a significant amount of fuctionality that FIELDS lacks, like functions to optimize the phasing arrangement of a cross section model, compute conductor heights needed to lower fields to target levels, and create a handful of different plots. For someone that knows a little Python, `emf.fields` is signifantly more flexible and more powerful than the old FIELDS program. More discussion of `emf.fields` can be found [here](http://mbaum1122.github.io/emf/README-emf.fields.html) and documentation of `emf.fields` can be found [here](http://mbaum1122.github.io/emf/emf.fields.html).

2. Like `emf.fields`, `emf.subcalc` was originally intended to supplement/extend the use of another modeling program. In this case, the other modeling program is [SUBCALC](http://www.enertech.net/html/emfw.html) (developed by [Enertech](http://www.enertech.net/html/emfw.html), sponsored by [EPRI](http://www.epri.com)). SUBCALC predicts magnetic fields over a fixed-height 2 dimensional grid and can model many non-parallel segments of power lines. It is generally used to model magnetic fields near electrical substations, where circuit arrangements can be complex. Just like `emf.fields`, `emf.subcalc` is now a complete replacement of the program it was originally intended to supplement, SUBCALC. Although it provides no fancy GUI, `emf.subcalc` is capable of calculating the magnetic fields produced by complex, three dimensional arrangements of current carrying wire segments. More discussion of `emf.subcalc` can be found [here](http://mbaum1122.github.io/emf/README-emf.subcalc.html) and documentation of `emf.subcalc` can be found [here](http://mbaum1122.github.io/emf/emf.subcalc.html).

###### Other Things

`emf` has mostly been used with Python 2.7. At different points it has been used compatibly with Python 3.x, but compatibility isn't routinely checked for.

###### Python Package Dependencies
* `os`
* `copy`
* `glob`
* `numpy`
* `scipy`
* `ctypes`
* `shutil`
* `pandas`
* `textwrap`
* `itertools`
* `matplotlib`
