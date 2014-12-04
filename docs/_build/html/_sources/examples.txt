examples
========

Quick Demo
----------

...

Calculate Blending Statistics for CFHT, DES, LSST
-------------------------------------------------

Run i-band calculations::

	./galsimcat.py -i OneDegSq.dat -x 0.5 -y 0.0 --max-size 30 --stamps --partials --save-field --save-noise --airmass 1.2 --extinction 0.07 -o lsst_i --pixel-scale 0.200 --width 4096 --height 4096 --exposure-time 6900 --sky-brightness 20.0 --zenith-fwhm 0.67 --zero-point 41.5 --hsm
	./galsimcat.py -i OneDegSq.dat -x 0.5 -y 0.0 --max-size 30 --stamps --partials --save-field --save-noise --airmass 1.2 --extinction 0.07 -o des_i  --pixel-scale 0.263 --width 3115 --height 3115 --exposure-time 1000 --sky-brightness 20.1 --zenith-fwhm 0.79 --zero-point 12.5 --hsm
	./galsimcat.py -i OneDegSq.dat -x 0.5 -y 0.0 --max-size 30 --stamps --partials --save-field --save-noise --airmass 1.2 --extinction 0.07 -o cfht_i --pixel-scale 0.185 --width 4428 --height 4428 --exposure-time 4300 --sky-brightness 20.3 --zenith-fwhm 0.64 --zero-point 10.0 --hsm

Run r-band calculations::

	./galsimcat.py -i OneDegSq.dat -x 0.5 -y 0.0 --max-size 30 --stamps --partials --save-field --save-noise --airmass 1.2 --extinction 0.10 -o lsst_r --pixel-scale 0.200 --width 4096 --height 4096 --exposure-time 6900 --sky-brightness 21.3 --zenith-fwhm 0.70 --zero-point 55.8 --hsm
	./galsimcat.py -i OneDegSq.dat -x 0.5 -y 0.0 --max-size 30 --stamps --partials --save-field --save-noise --airmass 1.2 --extinction 0.10 -o des_r  --pixel-scale 0.263 --width 3115 --height 3115 --exposure-time  800 --sky-brightness 21.1 --zenith-fwhm 0.79 --zero-point 16.8 --hsm
	./galsimcat.py -i OneDegSq.dat -x 0.5 -y 0.0 --max-size 30 --stamps --partials --save-field --save-noise --airmass 1.2 --extinction 0.10 -o cfht_r --pixel-scale 0.185 --width 4428 --height 4428 --exposure-time 2000 --sky-brightness 20.8 --zenith-fwhm 0.71 --zero-point 13.5 --hsm
