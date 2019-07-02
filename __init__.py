"""Weak lensing fast simulations and analysis for the LSST Dark Energy Science Collaboration.

This code was primarily developed to study the effects of overlapping sources on shear estimation,
photometric redshift algorithms, and deblending algorithms.
"""

__author__ = 'WeakLensingDeblending developers'
__email__ = 'dkirkby@uci.edu'
__version__ = '0.3dev'

import sys 
import os 

sys.path.insert(0, os.path.abspath(os.path.dirname(__file__)) )

import descwl 