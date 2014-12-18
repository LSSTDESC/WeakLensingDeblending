Simulation Output
=================

This document describes the content and format of the output file produced by the :ref:`prog-simulate` program.

The output consists of a `FITS file <http://fits.gsfc.nasa.gov/fits_primer.html>`_ containing several header/data units (HDUs).

You can inspect an output file's contents from an interactive python session, e.g.::

	from astropy.io import fits
	hdulist = fits.open('demo.fits')
	hdulist.info()
	hdulist[0].header
	hdulist.close()

Simulated Survey Image
----------------------

The primary HDU contains the final simulated image using double precision floats. All sources are superimposed in this source and fluxes are given in units of detected electrons during the full exposure time.

All of the :class:`descwl.survey.Survey` constructor args are saved as header keywords in the primary HDU, using only the first eight characters in upper case for the corresponding keys.

Analysis Results
----------------

HDU[1] contains a binary table where each row represents one simulated source.

You can load just the analysis results from the output file using, e.g.::

	import astropy.table
	table = astropy.table.Table.read('demo.fits',hdu=1)

To scroll through the table in an interactive python session, use::

	table.more()

To browse the table interactively (including seaching and sorting), use::

	table.show_in_browser(jsviewer=True)

To plot a histogram of signal-to-noise ratios for all visible galaxies (assuming that `matplotlib` is configured)::

	plt.hist(table['snr'][table['visible']])

Rendered Galaxy Stamps
----------------------

HDU[n+1] contains an image data cube for stamp n = 0,1,...  Each data cube HDU has header keywords `X_MIN` and `Y_MIN` that give the pixel offset of the stamp's lower-left corner from the lower-left corner of the full simulated survey image. Note that stamps may be partially outside of the survey image, but will always have some pixels above threshold within the image.

DS9 Usage
---------

If you open an output file with the `DS9 program <...>`_ you will normally only see the full simulated survey image in the primary HDU.  You can also use the `File > Open As > Multiple Extension Cube...` to view the nominal rendered stamp for each visible galaxy (but not any partial derivative images).
