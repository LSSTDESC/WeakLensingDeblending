"""Perform weak-lensing analysis of simulated sources.
"""

import numpy as np

import astropy.table

class OverlapAnalyzer(object):
    """Analyze impact of overlapping sources on weak lensing.

    Args:
        survey(descwl.survey.Survey): Simulated survey to describe with FITS header keywords.
    """
    def __init__(self,survey):
        self.survey = survey
        self.galaxies = [ ]

    def add_galaxy(self,model,stamps,bounds):
        """Add one galaxy to be analyzed.

        Args:
            model(:class:`descwl.model.Galaxy`): The galaxy model used for rendering.
            stamps(:class:`numpy.ndarray`): Array of shape (nstamp,width,height) containing
                pixel values for nstamp stamps of dimensions (width,height).
            bounds(galsim.BoundsI): Bounds of the stamps in the full simulated survey image.
        """
        self.galaxies.append((model,stamps,bounds))

    def finalize(self):
        """Finalize analysis of all added galaxies.

        Returns:
            :class:`astropy.table.Table`: Table of analysis results with one row per galaxy.
        """
        data = np.empty(len(self.galaxies),dtype=[
            ('visible',bool),('snr',float)
            ])
        for index,(model,stamps,bounds) in enumerate(self.galaxies):
            # Is this galaxy's centroid visible in the survey image?
            data['visible'][index] = self.survey.image.bounds.includes(bounds.center())
            # Calculate the SNR this galaxy would have without any overlaps.
            data['snr'][index] = 1.23

        return astropy.table.Table(data,copy = False)
