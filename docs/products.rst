Data Products
=============

The following data products related to this package are available for `download from SLAC <ftp://ftp.slac.stanford.edu/groups/desc/WL/>`_:

======================= ========== ========================================================
Filename                Size (Mb)  Description
======================= ========== ========================================================
OneDegSq.fits.gz              81.4 Input LSST DM galaxy catalog covering 1 sq.deg.
LSST_i.fits.gz               235.0 Simulated LSST i-band at full depth covering 1 chip
LSST_r.fits.gz               282.9 Simulated LSST r-band at full depth covering 1 chip
DES_i.fits.gz                 77.4 Simulated DES i-band at full depth covering same area
DES_r.fits.gz                 85.0 Simulated DES r-band at full depth covering same area
LSST_i_trimmed.fits.gz        27.3 Trimmed simulation results without per-object datacubes
LSST_r_trimmed.fits.gz        30.8 Trimmed simulation results without per-object datacubes
DES_i_trimmed.fits.gz         11.9 Trimmed simulation results without per-object datacubes
DES_r_trimmed.fits.gz         12.9 Trimmed simulation results without per-object datacubes
======================= ========== ========================================================

All files are in compressed FITS format and should be uncompressed after downloading::

	gunzip *.fits.gz

The commands used to generate the `OneDegSq.fits` catalog are :ref:`described here <catalog-create>`. The following commands were used to run the LSST and DES simulations::

	nohup ./simulate.py --catalog-name OneDegSq.fits --survey-name LSST --filter-band i --output-name LSST_i > LSST_i.log &
	nohup ./simulate.py --catalog-name OneDegSq.fits --survey-name LSST --filter-band r --output-name LSST_r > LSST_r.log &
	nohup ./simulate.py --catalog-name OneDegSq.fits --survey-name DES --filter-band i --output-name DES_i > DES_i.log &
	nohup ./simulate.py --catalog-name OneDegSq.fits --survey-name DES --filter-band r --output-name DES_r > DES_r.log &

The resulting FITS files were then trimmed to remove the individual object datacubes using::

	import astropy.io.fits
	for name in ('LSST_i','LSST_r','DES_i','DES_r'):
		hdus = astropy.io.fits.open(name+'.fits')
		del hdus[2:]
		hdus.writeto(name+'_trimmed.fits')
		hdus.close()
