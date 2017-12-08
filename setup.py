from setuptools import setup, find_packages

setup(
    name='emf',
    version='0.1',
    description='A program for modeling electric and magnetic fields (EMF) near parallel sets of power lines (2D models with emf.fields) and a program for modeling magnetic fields near complex arrangements of nonparellel power lines (3D models with emf.subcalc)',
    keywords='magnetic electric fields subcalc model modeling',
    author='Mark M. Baum',
    packages=['emf', 'emf.fields', 'emf.fields.enviro', 'emf.subcalc'],
    include_package_data=True,
    zip_safe=False,
    install_requires=['datetime',
            'numpy', 'pandas', 'matplotlib', 'scipy'],
            #'pkg_resources', 'itertools', 'textwrap', 'shutil', 'ctypes', 'copy', 'os'],
    classifiers=['Natural Language :: English',
            'Programming Language :: Python :: 3',
            'Topic :: Scientific/Engineering',
            'Topic :: Scientific/Engineering :: Physics']
)
