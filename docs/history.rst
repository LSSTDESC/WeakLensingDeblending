History
=======

The following changes were made between the original galsimcat.py and the new simulate.py:

- Fix catalog half-light radius bug: use sqrt(a*b) instead of catalog HLR value.
- Include AGN component when present in the catalog.
- Use Airy optical model instead of increasing atmospheric FWHM by 0.4 arcsec.
