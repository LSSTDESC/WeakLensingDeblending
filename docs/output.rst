Simulation Output
=================

This document describes the content and format of the output file produced by the :ref:`prog-simulate` program.

The output consists of a `FITS file <http://fits.gsfc.nasa.gov/fits_primer.html>`_ containing several header/data units (HDUs).

Primary HDU
-----------

The primary HDU[0] contains the final simulated image using double precision floats. All sources are superimposed in this source and fluxes are given in units of detected electrons during the full exposure time.

Extension HDUs
--------------

HDU[1] contains a binary table where each row represents one simulated source.

HDU[n+1] contains an image data cube for stamp n = 0,1,...
