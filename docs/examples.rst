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

Make a finder chart for overlapping groups::

	./display.py --input-name demo --annotate --annotate-format '%(grp_id)ld' --select 'grp_size>1' --select 'grp_rank==0' --magnification 2 --output-name finder.png

	./display.py --input-name demo --annotate --crop --group 402700184222 --magnification 16 --annotate-size x-large
	./fisher.py --input-name demo --galaxy 402700184222
	./fisher.py --input-name demo --group 402700184222
	./fisher.py --input-name demo --group 402700184222 --correlation

Calculate Blending Statistics for CFHT, DES, LSST
-------------------------------------------------

Run simulations to generate the :doc:`products` available for download from SLAC::

	./simulate.py --catalog-name OneDegSq.fits --survey-name LSST --filter-band i --output-name LSST_i
	./simulate.py --catalog-name OneDegSq.fits --survey-name LSST --filter-band r --output-name LSST_r
	./simulate.py --catalog-name OneDegSq.fits --survey-name DES --filter-band i --output-name DES_i
	./simulate.py --catalog-name OneDegSq.fits --survey-name DES --filter-band r --output-name DES_r

Each simulation takes about 30 minutes on a 2.5GHz i7 laptop and uses up to 5 Gb of memory.

Display a reduced simulated image::

	./display.py --input-name LSST_i --magnification 0.5 --output-name LSST_i.png
