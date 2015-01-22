To Do List
==========

* Formatted survey options printout.
* Use latest fitsio for writing simulation output.
* Use better 8-char names and add comments to FITS headers.
* Is the optical PSF making any difference? How much does it slow down simulate?
* Implement a 'exposure-time calculator' to calculate quantities for a profile specified on the command line (idea from Josh Meyers).
* Include a prior on dilation scale parameter for error estimation?
* Add ra,dec pointing and WCS to FITS primary HDU written by simulate.
* More flexible handling of the ra,dec catalog query window (we currently assume the query is centered near the origin).
* python 2-to-3 migration: from __future__ import absolute_import, division, print_function, unicode_literals
* Add zscale color gradient and angular scale legend options to display.py
* Add option to specify display view bounds independently of selected objects.
* Increase simulated area (from 1 LSST chip) for data products?

Longer-Term Projects
--------------------

* Validate galaxy input catalog quantities against existing datasets.
* Add stars from some catalog and estimate their contributions to blending statistics.
* Integrate over survey observing conditions (seeing, sky, ...) using operations simulator output.
* Test toy deblender performance on simulated images (Josh Myers?).
* Compare sextractor object detection with snr_grp detection threshold (Mandeep Gill).
* Compare DM pipeline object detection with snr_grp detection threshold.
* Use PhoSim to generate realistic bright-star templates that could be used in GalSim (Chris Walter?).
* Compare GalSim and PhoSim images generated for the same catalog footprint.
* Characterize galaxies that look like stars.
* Compare CFHTLS predictions with published results.
* Compare DES predictions with actual DES imaging.
* Add survey parameters for HSC to Survey class and compare predictions with actual HSC imaging (Rachel Mandelbaum?)
