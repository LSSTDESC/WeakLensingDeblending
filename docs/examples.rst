Examples
========

Command-line Options
--------------------

Print out usage info for command-line options::

	./simulate.py --help
	./display.py --help
	./fisher.py --help

Survey Parameters
-----------------

Print default camera and observing condition parameters for each supported (survey,filter) combination::

	./simulate.py --survey-defaults

Quick Simulation Demo
---------------------

Simulate an LSST i-band stack for a small (512x512) field::

	./simulate.py --catalog-name OneDegSq.fits --image-width 512 --image-height 512 --survey-name LSST --filter-band i --output-name demo --verbose

Make a finder chart for overlapping groups::

	./display.py --input-name demo --annotate --annotate-format '%(grp_id)ld' --select 'grp_size>1' --select 'grp_rank==0' --magnification 2 --output-name finder.png

Display a single blended group::

	./display.py --input-name demo --annotate --crop --group 402700184222 --magnification 16 --annotate-size x-large

Plot the Fisher matrix calculations for this group::

	./fisher.py --input-name demo --galaxy 402700184222
	./fisher.py --input-name demo --group 402700184222
	./fisher.py --input-name demo --group 402700184222 --correlation

Data Products Demo
------------------

Download the `LSST_i.fits` :doc:`data product <products>`, then display the location of galaxies that are unresolved due to blending::

	./display.py --input-name LSST_i --magnification 0.25 --select 'snr_grpf<5' --select 'snr_sky>10' --crosshair-color red
