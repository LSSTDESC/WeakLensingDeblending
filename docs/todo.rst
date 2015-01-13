To Do List
==========

* Formatted survey options printout.
* Use better 8-char names and add comments to FITS headers.
* Use transparency for selected objects in display.
* Add f/df,dg1,dg2 Fisher-matrix errors using isolated/blended assumption to output catalog.
* Is the optical PSF making any difference? How much does it slow down simulate?
* Implement a 'exposure-time calculator' to calculate quantities for a profile specified on the command line (idea from Josh).
* Include a prior on dilation scale parameter for error estimation?
* Add ra,dec pointing and WCS to FITS primary HDU written by simulate.
* More flexible handling of the ra,dec catalog query window (we currently assume the query is centered near the origin).
* Use print,div imports from future for python3.

Longer-Term Projects
--------------------

* Validate galaxy input catalog quantities against existing datasets.
* Add stars from some catalog and estimate their contributions to blending statistics.
* Integrate over survey observing conditions (seeing, sky, ...).
* Compare sextractor object detection with snr_grp detection threshold.
* Use PhoSim to generate realistic bright-star templates that could be used in GalSim.
* Compare GalSim and PhoSim images generated for the same catalog footprint.
* Add survey parameters for HSC.
* Compare CFHTLS predictions with published results.
