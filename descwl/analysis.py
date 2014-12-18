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

    def add_galaxy(self,stamps,x_min_pixels,y_min_pixels):
        """Add one galaxy to be analyzed.

        Args:
            stamps(:class:`numpy.ndarray`): Array of shape (nstamp,width,height) containing
                pixel values for nstamp stamps of dimensions (width,height).
            x_min_pixels(int): Left edge of stamps in pixels relative to the left edge
                of the simulated survey image.
            y_min_pixels(int): Bottom edge of stamps in pixels relative to the bottom edge
                of the simulated survey image.
        """
        self.stamps.append(stamps)
        nstamps,width,height = stamps.shape
        bounds = galsim.BoundsI(x_min_pixels,y_min_pixels,
            x_min_pixels+width-1,y_min_pixels+height-1)
        self.bounds.append(bounds)

    def finalize(self):
        """Finalize analysis of all added galaxies.

        Returns:
            :class:`astropy.table.Table`: Table of analysis results with one row per galaxy.
        """
        return None
