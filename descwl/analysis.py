"""Perform weak-lensing analysis of simulated sources.
"""

import galsim

class OverlapAnalyzer(object):
    """Analyze impact of overlapping sources on weak lensing.

    Args:
        survey(descwl.survey.Survey): Simulated survey to describe with FITS header keywords.
    """
    def __init__(self,survey):
        self.survey = survey
        self.stamps = [ ]
        self.bounds = [ ]

    def add_galaxy(self,stamps,bounds):
        """Add one galaxy to be analyzed.

        Args:
            stamps(:class:`numpy.ndarray`): Array of shape (nstamp,width,height) containing
                pixel values for nstamp stamps of dimensions (width,height).
            bounds(galsim.BoundsI): Bounds of the stamps in the full simulated survey image.
        """
        self.stamps.append(stamps)
        self.bounds.append(bounds)

    def finalize(self):
        """Finalize analysis of all added galaxies.

        Returns:
            :class:`astropy.table.Table`: Table of analysis results with one row per galaxy.
        """
        return None
