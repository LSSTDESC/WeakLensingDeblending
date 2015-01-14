Data Products
=============

The following data products related to this package are available for download from SLAC:

================= ========== =====================================================
Filename          Size (Mb)  Description
================= ========== =====================================================
OneDegSq.fits          130.5 Input LSST DM galaxy catalog covering 1 sq.deg.
LSST_i.fits             ~700 Simulated LSST i-band at full depth covering 1 chip
LSST_r.fits             ~700 Simulated LSST r-band at full depth covering 1 chip
DES_i.fits              ~500 Simulated DES i-band at full depth covering same area
DES_r.fits              ~500 Simulated DES r-band at full depth covering same area
================= ========== =====================================================

Details on how to download these files will be available soon.

The commands used to generate the `OneDegSq.fits` catalog are :ref:`described here <catalog-create>`. The following commands were used to run the LSST and DES simulations:

	nohup ./simulate.py --catalog-name OneDegSq.fits --survey-name LSST --filter-band i --output-name LSST_i > LSST_i.log &
	nohup ./simulate.py --catalog-name OneDegSq.fits --survey-name LSST --filter-band r --output-name LSST_r > LSST_r.log &
	nohup ./simulate.py --catalog-name OneDegSq.fits --survey-name DES --filter-band i --output-name DES_i > DES_i.log &
	nohup ./simulate.py --catalog-name OneDegSq.fits --survey-name DES --filter-band r --output-name DES_r > DES_r.log &

The result FITS files were then trimmed to remove individual object datacubes using::

	import astropy.io.fits
	for name in ('LSST_i','LSST_r','DES_i','DES_r'):
		hdus = astropy.io.fits.open(name+'.fits')
		del hdus[2:]
		hdus.writeto(name+'_trimmed.fits')
		hdus.close()
