from setuptools import setup

setup(
    name='emf',
    version='0.1',
    description='A program for modeling electromagnetic fields (EMF) near parallel sets of power lines (emf.fields) and a program for analyzing the results of another EMF model called SubCalc (emf.subcalc)',
    url='https://github.com/mbaum1122/emf',
    download_url='https://github.com/mbaum1122/emf',
    author='Mark M. Baum',
    author_email='mbaum1122@gmail.com',
    maintainer='Mark M. Baum',
    maintainer_email='mbaum1122@gmail.com',
    packages=['emf', 'emf.fields', 'emf.subcalc'],
    include_package_data=True,
    zip_safe=False,
    license='MIT',
    install_requires=['os', 'copy', 'glob', 'shutil', 'itertools', 'numpy',
                            'pandas', 'matplotlib', 'scipy', 'ctypes'],
    classifiers=['Natural Language :: English',
                'Programming Language :: Python :: 2.7',
                'Topic :: Scientific/Engineering',
                'Topic :: Scientific/Engineering :: Physics']
)
