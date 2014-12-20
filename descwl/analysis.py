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
            ('ab_mag',np.float32),
            ('snr_sky',np.float32),
            ('snr_iso',np.float32),
            ])

        # Calculate isolated galaxy quantities and identify overlapping groups.
        data['grp_id'] = np.arange(num_galaxies)
        for index,(model,stamps,bounds) in enumerate(self.galaxies):
            data['db_id'][index] = model.identifier
            data['z'][index] = model.redshift
            data['ab_mag'][index] = model.ab_magnitude
            # Is this galaxy's centroid visible in the survey image?
            data['visible'][index] = self.survey.image.bounds.includes(bounds.center())
            # Calculate the SNR this galaxy would have without any overlaps in the
            # sky-dominated limit.
            sky = self.survey.mean_sky_level
            fiducial = stamps[0].flatten()
            data['snr_sky'][index] = np.sqrt(np.sum(fiducial**2)/sky)
            # Include the signal in the noise variance.
            mu0 = fiducial + sky
            data['snr_iso'][index] = np.sqrt(np.sum(fiducial**2*(mu0**-1 + 0.5*mu0**-2)))
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
        num_groups = np.max(data['grp_id']) + 1
        for grp_id in range(num_groups):
            grp_members = (data['grp_id'] == grp_id)
            grp_size = np.count_nonzero(grp_members)
            data['grp_size'][grp_members] = grp_size

        return astropy.table.Table(data,copy = False)