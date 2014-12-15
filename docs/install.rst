Installation
============

Install with GIT
----------------

The code is hosted on `github <https://github.com/DarkEnergyScienceCollaboration/WeakLensingDeblending>`_ so the easiest method to perform an initial installation is with `git <http://git-scm.com>`_::

	git clone git@github.com:DarkEnergyScienceCollaboration/WeakLensingDeblending.git

This will create a new subdirectory `WeakLensingDeblending` containing the latest stable version of the complete package.

Update with GIT
---------------

You can update your local copy of the package at any time using::

	cd WeakLensingDeblending
	git update

Getting Started
---------------

Programs can be run directly from the top-level directory without needing to set `PYTHONPATH` as long as you have the required packages already installed, e.g.::

	cd WeakLensingDeblending
	./simulate.py --help

For an introduction to the available programs, see :doc:`here </programs>` and for examples of running these programs see :doc:`here </examples>`.

Required Packages
-----------------

The following python packages are required by this package:

* numpy
* astropy
* `galsim <https://github.com/GalSim-developers/GalSim>`_

Note that, except for galsim, these packages are all available in recent `anaconda <https://store.continuum.io/cshop/anaconda/>`_ or `enthought canopy <https://www.enthought.com/products/canopy/>`_ distributions.
