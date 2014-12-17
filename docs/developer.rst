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
