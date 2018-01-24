"""Weak lensing fast simulations and analysis for the LSST Dark Energy Science Collaboration.

This code was primarily developed to study the effects of overlapping sources on shear estimation,
photometric redshift algorithms, and deblending algorithms.
"""

__author__ = 'WeakLensingDeblending developers'
__email__ = 'dkirkby@uci.edu'
__version__ = '0.3dev'

from . import catalog
from . import survey
from . import model
from . import render
from . import analysis
from . import output
from . import trace
