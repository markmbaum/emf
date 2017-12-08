# `emf`

![equations](docs/img/both-equations.png)

The `emf` package is a container for two subpackages (documentation accessible [here](http://mbaum1122.github.io/emf/)). The subpackages, called `fields` and `subcalc`, are for modeling electric and magnetic fields near power lines. The `fields` subpackage is for 2D models. It predicts electric and magnetic fields near sets of parallel power lines. The `subcalc` package is for 3D models, predicting the magnetic fields generated by complex, nonparallel collections of current carrying wire segments. Both subpackages can be used to build, compute, and analyze models without any other programs. They rely heavily on `numpy`, `pandas`, and `matplotlib`, but the EMF calculations do not rely on algorithms from any external libraries or legacy code. For 2D models in `emf.fields`, the fields calculations are done only in Python. For 3D models in `emf.subcalc`, a [C extension](emf/subcalc/lift/lift.c) is used to speed up the intensive number crunching.

Both `fields` and `subcalc` were originally meant to supplement older modeling programs but eventually became complete replacements for them. Thus, both packages are styled after the older programs and contain a fair amount of functionality for working with those programs. Both `fields` and `subcalc` also have specially formatted excel/csv templates for storing input data. Although these template files are never strictly necessary, they make the packages much more usable for those who only know a little bit of Python and they make larger-scale modeling projects more repeatable and less error-prone. Neither package has a UI, just scripts and (optionally) template files.

### `emf.fields` resources

* [Further discussion](docs/README-fields.md)
* [Documentation](http://mbaum1122.github.io/emf/emf.fields.html)
* [Notebook - `fields` modeling from scratch](docs/notebooks/fields/fields-workflow-from-scratch.ipynb)
* [Notebook - optimizing underground circuits](docs/notebooks/fields/underground-line-optimization.ipynb)
* [Notebook - using a `fields` template](docs/notebooks/fields/using-a-template.ipynb)
* [Notebook - varying the depth of an underground circuit](docs/notebooks/fields/underground-delta-depth-test.ipynb)

### `emf.subcalc` resources

* [Further discussion](docs/README-subcalc.md)
* [Documentation](http://mbaum1122.github.io/emf/emf.subcalc.html)
* [Notebook - build a model, compute results, and make plots](docs/notebooks/subcalc/small-model-tutorial.ipynb)
* [Notebook - using `emf.subcalc` templates](docs/notebooks/subcalc/tower-and-footprint-templates.ipynb)
* [Notebook - sampling `Model` objects](docs/notebooks/subcalc/sampling-model-objects.ipynb)
* [Notebook - working with SUBCALC input and output files](docs/notebooks/subcalc/working-with-SUBCALC-files.ipynb)

### Installation

To use `emf` you must first install Python 3 and `emf`'s dependencies. `emf` relies on many other widely available Python libraries, which are listed below. [Anaconda](https://www.continuum.io/downloads) is a widely used distribution that should include everything you need to use `emf`.

Once Python and the relevant libraries are installed, you can install `emf` by running the `setup.py` script in the 'emf' package's directory. In a command prompt or terminal, simply navigate to the directory that contains `setup.py` and enter the following command.

```python
python setup.py install
```

If you download an updated version of `emf` or otherwise modify the package, you'll need to reinstall it by rerunning `setup.py` with the same command as above.

##### Python Package Dependencies
`os`, `copy`, `numpy`, `scipy`, `ctypes`, `shutil`, `pandas`, `datetime`, `textwrap`, `itertools`, `matplotlib`, `pkg_resources`
