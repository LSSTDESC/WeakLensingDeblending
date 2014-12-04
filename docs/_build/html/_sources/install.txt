install
=======

Programs can be run directly from the top-level directory without needing to set `PYTHONPATH` as long as you have the required packages already installed, e.g.::

	git clone git@github.com:DarkEnergyScienceCollaboration/WeakLensingDeblending.git
	cd WeakLensingDeblending
	./galsimcat.py --help

Required Packages
-----------------

The following python packages are required by this package:

* numpy
* `galsim <https://github.com/GalSim-developers/GalSim>`_
