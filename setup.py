#!/usr/bin/env python
# -*- coding: utf-8 -*-

from setuptools import setup

requirements = [
    'fitsio',
    'galsim',
    'numpy',
    'astropy',
    'lmfit',
    'six'
]

setup(
    name='descwl',
    version='0.3dev',
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
