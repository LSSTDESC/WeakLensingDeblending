Revision History
================

v0.?
----
- Version used for DESC highlight project data products.
- Refactor code and add sphinx documentation.
- Rename `galsimcat.py` as `simulate.py` and `lsst2wl.py` to `dbquery.py`.
- Add new programs `display.py`, `fisher.py`.
- Fix catalog half-light radius bug: use sqrt(a*b) instead of catalog HLR value.
- Include AGN component when present in the catalog.
- Use Airy optical model instead of increasing atmospheric FWHM by 0.4 arcsec.
- Implement new blending figures of merit.

v0.1
----
- Version used to prepare results for the `Dec 2013 LSST DESC meeting <https://indico.bnl.gov/conferenceDisplay.py?confId=691>`_ in Pittsburgh.
