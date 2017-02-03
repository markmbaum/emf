# `emf`

![equations](docs/img/both-equations.png)

The `emf` package is a container for two subpackages (documentation accessible [here](http://mbaum1122.github.io/emf/)). The subpackages, which are called `fields` and `subcalc`, are for modeling electric and magnetic fields near power lines. The `fields` subpackage is for 2D models. It predicts electric and magnetic fields near sets of parallel power lines. The `subcalc` package is for 3D models, predicting the magnetic fields generated by complex, nonparallel collections of current carrying wire segments. Both subpackages are freestanding, meaning they can be used to build, compute, and analyze models all by themselves without any other programs. They rely heavily on `numpy`, `pandas`, and `matplotlib`, but the EMF calculations do not rely on algorithms from any external libraries or legacy code.

Both `fields` and `subcalc` were originally meant to supplement older modeling programs and eventually became complete replacements of them. Thus, both packages are styled after the older programs and still contain a fair amount of functionality for working with those programs. Both packages also have specially formatted excel/csv files for storing input data. Although these template files are never strictly necessary, they make the packages much more usable for those who only know a little bit of Python and they make larger-scale modeling projects more repeatable and better documented. Neither package has a fancy UI, just scripts and (optionally) template files.

### `emf.fields` resources

* [Further discussion](http://mbaum1122.github.io/emf/README-fields.html)
* [Documentation](http://mbaum1122.github.io/emf/emf.fields.html)
* [Tutorial - `fields` Modeling From Scratch](https://github.com/mbaum1122/emf/blob/master/docs/notebooks/fields/fields-workflow-from-scratch.ipynb)
* [Example - Optimizing Underground Circuits](https://github.com/mbaum1122/emf/blob/master/docs/notebooks/fields/underground-line-optimization.ipynb)
* [Example - Using a `fields` Template](https://github.com/mbaum1122/emf/blob/master/docs/notebooks/fields/using-a-template.ipynb)

### `emf.subcalc` resources

* [Further discussion](http://mbaum1122.github.io/emf/README-subcalc.html)
* [Documentation](http://mbaum1122.github.io/emf/emf.subcalc.html)
* [Tutorial - build a model, compute results, and make plots](https://github.com/mbaum1122/emf/blob/master/docs/notebooks/subcalc/small-model-tutorial.ipynb)
* [Example - using `emf.subcalc` templates](https://github.com/mbaum1122/emf/blob/master/docs/notebooks/subcalc/tower-and-footprint-templates.ipynb)

##### Other Things

`emf` has mostly been used with Python 2.7. At different points it has been used compatibly with Python 3.x, but compatibility hasn't been checked for a while.

##### Python Package Dependencies
`os`, `copy`, `glob`, `numpy`, `scipy`, `ctypes`, `shutil`, `pandas`, `textwrap`, `itertools`, `matplotlib`
