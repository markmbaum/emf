from setuptools import setup

setup(
    name='emf',
    version='0.1',
    description='A program for modeling electric and magnetic fields (EMF) near parallel sets of power lines (2D models with emf.fields) and a program for modeling magnetic fields near complex arrangements of nonparellel power lines (3D models with emf.subcalc)',
    keywords='magnetic electric fields subcalc model modeling',
    author='Mark M. Baum',
    packages=['emf', 'emf.fields', 'emf.subcalc'],
    include_package_data=True,
    zip_safe=False,
    install_requires=['os', 'copy', 'ctypes', 'shutil', 'datetime', 'textwrap',
            'itertools', 'numpy', 'pandas', 'matplotlib', 'scipy',
            'pkg_resources'],
    classifiers=['Natural Language :: English',
            'Programming Language :: Python :: 2.7',
            'Topic :: Scientific/Engineering',
            'Topic :: Scientific/Engineering :: Physics']
)
