#!/usr/bin/env python
import distutils
from distutils.core import setup

description = "Fast simulations and analysis for LSST DESC"

setup(name="descwl", 
      version="0.1",
      description=description,
      url="https://github.com/LSSTDESC/WeakLensingDeblending",
      author="David Kirkby, Javier Sanchez, Ismael Mendoza LSST DESC",
      author_email="dkirkby@uci.edu, francs1@uci.edu",
      packages=['descwl'])
