.. WeakLensingDeblending documentation master file, created by
   sphinx-quickstart on Wed Dec  3 17:14:11 2014.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. image:: img/descwl.*

DESC Weak Lensing Deblending Package
====================================

Fast simulations and analysis for the Weak Lensing Working Group of the LSST `Dark Energy Science Collaboration <http://www.lsst-desc.org>`_.

This software was primarily developed to study the effects of overlapping sources on shear estimation,
photometric redshift algorithms, and deblending algorithms. Users can run their own simulations (of LSST and other surveys) or, else download the :doc:`galaxy catalog<catalog>` and :doc:`simulation outputs<products>` to use with their own code or analyze with the tools provided here.

The code is hosted on `github <https://github.com/DarkEnergyScienceCollaboration/WeakLensingDeblending>`_.  Please use the `issue tracker <https://github.com/DarkEnergyScienceCollaboration/WeakLensingDeblending/issues>`_ to let us know about any issues you have with installing or running this code, or to request new features.

.. toctree::
   :maxdepth: 2

   install
   catalog
   examples
   programs
   output
   products
   notebooks
   todo
   developer
   history

Modules API Reference
---------------------

.. toctree::
   :maxdepth: 3

   src/descwl
