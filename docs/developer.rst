Information for Developers
==========================

Build the documentation
-----------------------

To build a local copy of the HTML documentation, use::

	cd .../descwl/docs
	make html

To view the most recent build of the HTML documentation, point your browser to `.../docs/_build/html/index.html`

To create a tarball snapshot `.../descwl/docs/_build/descwl.tgz` that can be installed on a web server, use::

	cd .../descwl/docs/_build/html
	tar -zcf ../descwl.tgz .

Add a new package module
------------------------

Create a new file `descwl/xxx.py` for module `descwl.xxx` with an initial descriptive docstring.

Add the line `import xxx` to `descwl/__init__.py`.

Add the line `descwl.xxx` to the list of submodules in `docs/src/descwl.rst`.

Create a new documentation file `docs/src/descwl.xxx.rst` containing (replace `xxx` in two places)::

	descwl.xxx module
	=================

	.. automodule:: descwl.xxx
	    :members:
	    :undoc-members:
	    :show-inheritance:

Profiling
---------

The `simulate` program has a `--memory-trace` option that displays the memory usage at frequent checkpoints during the execution.

Profile CPU usage with the `standard recipe <https://docs.python.org/2/library/profile.html#instant-user-s-manual>`_, for example::

	python -m cProfile -o profile.out ./simulate.py --catalog-name OneDegSq.fits --output-name profile

Examine the results using, for example::

	import pstats
	p = pstats.Stats('profile.out')
	p.sort_stats('time').print_stats(10)

Here is a sample profiling output using the above examples::

	Sat Jan  3 12:16:55 2015    profile.out

	         217334792 function calls (216127873 primitive calls) in 813.478 seconds

	   Ordered by: internal time
	   List reduced from 2518 to 10 due to restriction <10>

	   ncalls  tottime  percall  cumtime  percall filename:lineno(function)
	    50963  520.836    0.010  569.392    0.011 /Users/david/anaconda/lib/python2.7/site-packages/galsim/base.py:873(drawImage)
	        1   41.916   41.916   48.871   48.871 ./descwl/analysis.py:142(finalize)
	   593789   40.541    0.000   42.049    0.000 /Users/david/anaconda/lib/python2.7/site-packages/galsim/base.py:212(__init__)
	   281014   15.681    0.000   15.681    0.000 {method 'reduce' of 'numpy.ufunc' objects}
	    51547   15.326    0.000  611.955    0.012 ./descwl/render.py:71(render_galaxy)
	  2946402    9.266    0.000   37.495    0.000 /Users/david/anaconda/lib/python2.7/site-packages/astropy/config/configuration.py:375(__call__)
	 41011977    8.557    0.000    8.977    0.000 {isinstance}
	  3457423    8.052    0.000   13.553    0.000 /Users/david/anaconda/lib/python2.7/site-packages/astropy/config/configuration.py:622(get_config)
	  1178324    4.791    0.000   21.044    0.000 /Users/david/anaconda/lib/python2.7/site-packages/astropy/io/fits/card.py:555(value)
	   367098    4.396    0.000    5.389    0.000 /Users/david/anaconda/lib/python2.7/site-packages/galsim/image.py:188(__init__)
