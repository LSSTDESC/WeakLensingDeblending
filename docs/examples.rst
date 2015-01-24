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

	./display.py --input-name demo --info '%(grp_id)ld' --select 'grp_size>1' --select 'grp_rank==0' --magnification 2 --output-name finder.png

Display a single blended group::

	./display.py --input-name demo --info 'z=%(z).1f' --crop --group 2202242785 --magnification 16

Plot the Fisher matrix calculations for this group::

	./fisher.py --input-name demo --galaxy 2202242785 --verbose
	./fisher.py --input-name demo --group 2202242785 --verbose
	./fisher.py --input-name demo --group 2202242785 --correlation

Data Products Demo
------------------

Download the `LSST_i.fits` :doc:`data product <products>`, then display the location of galaxies that are unresolved due to blending::

	./display.py --input-name LSST_i --magnification 0.25 --select 'snr_grpf<5' --select 'snr_sky>10' --crosshair-color red

Compare with SExtractor detected objects for one cluster::

	./display.py --input-name LSST_i --crop --group 402700079754 --magnification 8 --match-catalog LSST_i_se.cat

Compare simulations of the same region in different surveys::

	./display.py -i LSST_i --magnification 10 --select-region '[-10.20,2.80,-30.40,-12.80]' --view-region '[-10.20,2.80,-30.40,-12.80]' --info '%(ab_mag).1f\n%(z).1f' --info-color black --info-size x-large --hide-background -o LSST_cat.png

	./display.py -i LSST_i --match-catalog LSST_i_noise.cat --magnification 10 --select-region '[-10.20,2.80,-30.40,-12.80]' --view-region '[-10.20,2.80,-30.40,-12.80]' --add-noise 1 --info '%(purity).2f\n%(snr_grpf).1f(%(snr_sky).1f)' --info-color red --info-size x-large --match-color yellow --crosshair-color red --hide-selected -o LSST_sim.png

	./display.py -i DES_i --match-catalog DES_i_noise.cat --magnification 13.15 --select-region '[-10.20,2.80,-30.40,-12.80]' --view-region '[-10.20,2.80,-30.40,-12.80]' --add-noise 1 --info '%(purity).2f\n%(snr_grpf).1f(%(snr_sky).1f)' --info-color red --info-size x-large --match-color yellow --crosshair-color red --hide-selected -o DES_sim.png
