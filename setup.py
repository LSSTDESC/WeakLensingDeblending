#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
from setuptools import setup

requirements = [
    'fitsio',
    'galsim',
    'numpy',
    'astropy',
    'lmfit',
    'six'
]

pth = os.path.join(
    os.path.dirname(os.path.realpath(__file__)),
    "descwl",
    "version.py"
)
with open(pth, 'r') as fp:
    exec(fp.read())

setup(
    name='descwl',
    version=__version__,
    description='Weak lensing fast simulations and analysis for the LSST DESC',
    long_description='Weak lensing fast simulations and analysis for the LSST DESC',
    author='descwl developers',
    author_email='dkirkby@uci.edu',
    url='https://github.com/LSSTDESC/WeakLensingDeblending',
    packages=[
        'descwl',
    ],
    package_dir={'descwl': 'descwl'},
    scripts = [ ],
    #include_package_data=True,
    #zip_safe=False,
    install_requires=[ ], #requirements,
    license='MIT',
)
