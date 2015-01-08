Examples
========

Command-line Options
--------------------

Print out usage info for command-line options::

	./simulate.py --help

Survey Parameters
-----------------

Print default camera and observing condition parameters for each supported (survey,filter) combination::

	./simulate.py --survey-defaults

Quick Demo
----------

Run i-band calculation for LSST with a small field and verbose output::

	./simulate.py --catalog-name OneDegSq.fits --image-width 512 --image-height 512 --survey-name LSST --filter-band i --output-name demo --verbose --verbose-model --verbose-render

Make a finder chart for overlapping groups:

	./display.py --input-name demo --annotate --annotate-format '%(grp_id)ld' --select 'grp_size>1' --select 'grp_rank==0' --magnification 2 --output-name finder.png

	./display.py --input-name demo --annotate --crop --group 102363 --magnification 16 --verbose
	./fisher.py --input-name demo --galaxy 102363
	./fisher.py --input-name demo --group 102363
	./fisher.py --input-name demo --group 102363 --correlation

Calculate Blending Statistics for CFHT, DES, LSST
-------------------------------------------------

Run i-band calculations (each one takes about 30 minutes on a 2.5GHz i7 laptop and uses up to 5 Gb of memory)::

	./simulate.py --catalog-name OneDegSq.fits --survey-name LSST --filter-band i --output-name LSST_i
	./simulate.py --catalog-name OneDegSq.fits --survey-name DES --filter-band i --output-name DES_i

	./display.py --input-name LSST_i --magnification 0.5 --output-name LSST_i.png

	./galsimcat.py -i OneDegSq.dat -x 0.5 -y 0.0 --max-size 30 --stamps --partials --save-field --save-noise --airmass 1.2 --extinction 0.07 -o lsst_i --pixel-scale 0.200 --width 4096 --height 4096 --exposure-time 6900 --sky-brightness 20.0 --zenith-fwhm 0.67 --zero-point 41.5 --hsm
	./galsimcat.py -i OneDegSq.dat -x 0.5 -y 0.0 --max-size 30 --stamps --partials --save-field --save-noise --airmass 1.2 --extinction 0.07 -o des_i  --pixel-scale 0.263 --width 3115 --height 3115 --exposure-time 1000 --sky-brightness 20.1 --zenith-fwhm 0.79 --zero-point 12.5 --hsm
	./galsimcat.py -i OneDegSq.dat -x 0.5 -y 0.0 --max-size 30 --stamps --partials --save-field --save-noise --airmass 1.2 --extinction 0.07 -o cfht_i --pixel-scale 0.185 --width 4428 --height 4428 --exposure-time 4300 --sky-brightness 20.3 --zenith-fwhm 0.64 --zero-point 10.0 --hsm

Run r-band calculations::

	./galsimcat.py -i OneDegSq.dat -x 0.5 -y 0.0 --max-size 30 --stamps --partials --save-field --save-noise --airmass 1.2 --extinction 0.10 -o lsst_r --pixel-scale 0.200 --width 4096 --height 4096 --exposure-time 6900 --sky-brightness 21.3 --zenith-fwhm 0.70 --zero-point 55.8 --hsm
	./galsimcat.py -i OneDegSq.dat -x 0.5 -y 0.0 --max-size 30 --stamps --partials --save-field --save-noise --airmass 1.2 --extinction 0.10 -o des_r  --pixel-scale 0.263 --width 3115 --height 3115 --exposure-time  800 --sky-brightness 21.1 --zenith-fwhm 0.79 --zero-point 16.8 --hsm
	./galsimcat.py -i OneDegSq.dat -x 0.5 -y 0.0 --max-size 30 --stamps --partials --save-field --save-noise --airmass 1.2 --extinction 0.10 -o cfht_r --pixel-scale 0.185 --width 4428 --height 4428 --exposure-time 2000 --sky-brightness 20.8 --zenith-fwhm 0.71 --zero-point 13.5 --hsm
