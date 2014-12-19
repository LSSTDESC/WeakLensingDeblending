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
        # Define columns and allocate space for our table data.
        num_galaxies = len(self.galaxies)
        data = np.empty(num_galaxies,dtype=[
            ('db_id',np.int64),
            ('grp_id',np.int16),
            ('grp_size',np.int16),
            ('visible',np.bool8),
            ('z',np.float32),
            ('snr_sky',np.float32)
            ])

        data['grp_id'] = np.arange(num_galaxies)
        for index,(model,stamps,bounds) in enumerate(self.galaxies):
            data['db_id'][index] = model.identifier
            data['z'][index] = model.redshift
            # Is this galaxy's centroid visible in the survey image?
            data['visible'][index] = self.survey.image.bounds.includes(bounds.center())
            # Calculate the SNR this galaxy would have without any overlaps in the
            # sky-dominated limit.
            flat_fiducial = stamps[0].flat
            data['snr_sky'][index] = np.sqrt(
                np.dot(flat_fiducial,flat_fiducial)/self.survey.mean_sky_level)
            # Loop over earlier galaxies to build overlapping groups.
            for pre_index,(pre_model,pre_stamps,pre_bounds) in enumerate(self.galaxies[:index]):
                # Do bounding boxes overlap?
                overlap = bounds & pre_bounds
                if overlap.area() == 0:
                    continue
                # Are there any overlapping pixels with non-zero flux from both sources?
                overlap_flux_product = np.sum(
                    stamps[0,
                        overlap.ymin-bounds.ymin:overlap.ymax-bounds.ymin+1,
                        overlap.xmin-bounds.xmin:overlap.xmax-bounds.xmin+1]*
                    pre_stamps[0,
                        overlap.ymin-pre_bounds.ymin:overlap.ymax-pre_bounds.ymin+1,
                        overlap.xmin-pre_bounds.xmin:overlap.xmax-pre_bounds.xmin+1])
                if overlap_flux_product == 0:
                    continue
                # Decide which group joins the other.
                grp_id_new = min(data['grp_id'][index],data['grp_id'][pre_index])
                grp_id_old = max(data['grp_id'][index],data['grp_id'][pre_index])
                # Reassign all galaxies in grp_id_old to grp_id_new.
                data['grp_id'][data['grp_id'] == grp_id_old] = grp_id_new

        # Analyze overlapping groups.
        data['grp_size'] = 0

        return astropy.table.Table(data,copy = False)
