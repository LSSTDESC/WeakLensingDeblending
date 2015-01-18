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

* numpy (version >= 1.9)
* astropy (version >= 0.4)
* `fitsio <https://github.com/esheldon/fitsio>`_ (version >= 0.9.6)
* `galsim <https://github.com/GalSim-developers/GalSim>`_ (version >= 1.2)

Note that `numpy` and `astropy` are both available in recent `anaconda <https://store.continuum.io/cshop/anaconda/>`_ or `enthought canopy <https://www.enthought.com/products/canopy/>`_ distributions. The `fitsio` package is required for performance reasons, although it overlaps the functionality of the pure-python `astropy.io.fits` module. Installing GalSim is a more involved process, but well worth the effort.

You can check your astropy version using::

	import astropy.version
	print astropy.version.version

A version of at least 0.4 is required due to recent changes in tables and FITS I/O. If you have a pre-0.4 version of astropy via anaconda, you can update using::

	sudo conda update conda
	sudo conda update anaconda

Note that some additional packages are required to :ref:`query the LSST DM catalog <catalog-create>`, but you will not normally need to do this.

The `psutil <https://pypi.python.org/pypi/psutil>`_ package is required if you use the `--memory-trace` command-line argument to the :ref:`prog-simulate` program, but you normally would not need to do this.
