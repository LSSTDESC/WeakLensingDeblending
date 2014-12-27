"""Perform weak-lensing analysis of simulated sources.
"""

import numpy as np

import astropy.table

import galsim

class OverlapResults(object):
    """Results of analyzing effects of overlapping sources on weak lensing.

    Results are normally obtained by calling :meth:`OverlapAnalyzer.finalize` after
    simulating an image or using :meth:`descwl.output.Reader` to load saved
    simulation results.

    Indices returned by the `find` methods below can be used to lookup analysis
    results via `table[index]` and the corresponding simulated datacube using
    `stamps[index]` and `bounds[index]`.

    Args:
        survey(descwl.survey.Survey): Simulated survey that results are based on.
        table(astropy.table.Table): Table of analysis results with one row per galaxy.
        stamps(array): Array of :class:`numpy.ndarray` postage-stamp datacubes.
        bounds(array): Array of galsim.BoundsI objects giving pixel coordinates of
            each datacube in the full simulated image.
    """
    def __init__(self,survey,table,stamps,bounds):
        self.survey = survey
        self.table = table
        self.stamps = stamps
        self.bounds = bounds
        self.num_objects = len(self.table)
        self.locals = { name: self.table[name] for name in self.table.colnames }

    def select(self,selector):
        """Select objects.

        This function is implemented using :func:`eval` so should only be used when the
        source of the selector string is trusted.

        Args:
            selector(str): A selection string consisting of a boolean expression of
                column table names, e.g. 'snr_iso > 5'.  The special strings 'ALL' and
                'NONE' do what you would expect.

        Returns:
            :class:`numpy.ndarray`: Boolean array containing a selection mask with one
                entry per object whose True/False value indicates if the object has
                been selected.

        Raises:
            RuntimeError: Selector uses an undefined variable name.
        """
        if selector == 'ALL':
            return np.ones(self.num_objects,dtype=bool)
        elif selector == 'NONE':
            return np.zeros(self.num_objects,dtype=bool)
        else:
            try:
                return eval(selector,self.locals)
            except NameError,e:
                raise RuntimeError('%s in selector %r.' % (e.message,selector))

    def get_stamp(self,index,datacube_index=0):
        """Return the simulated postage stamp for a single galaxy.

        Args:
            index(int): Index of the requested galaxy in our results. Use
                meth:`find_galaxy` to get the index for an identifier.
            datacube_index(int): Which slice of the datacube to return.

        Returns:
            galsim.ImageView: A GalSim image for the requested galaxy.

        Raises:
            RuntimeError: No such galaxy in the results.
        """
        if index < 0 or index >= len(self.stamps):
            raise RuntimeError('No such galaxy with index=%d.' % index)
        bounds = self.bounds[index]
        image = galsim.Image(self.stamps[index][datacube_index],
            xmin = bounds.xmin,ymin = bounds.ymin,
            scale = self.survey.pixel_scale,make_const = True)
        return image

    def get_subimage(self,indices):
        """Return simulated subimage of a set of objects.

        Args:
            indices(iterable): Indices of the objects to include in the subimage.

        Returns:
            galsim.Image: Image of the selected objects or None of indices is empty.

        Raises:
            RuntimeError: An index is out of range.
        """
        if len(indices) == 0:
            return None
        subimage_bounds = galsim.BoundsI()
        try:
            for index in indices:
                subimage_bounds += self.bounds[index]
        except IndexError:
            raise RuntimeError('Index out of range %d' % index)
        subimage = galsim.Image(bounds = subimage_bounds, scale = self.survey.pixel_scale)
        for index in indices:
            stamp = self.get_stamp(index)
            overlap = subimage_bounds & stamp.bounds
            subimage[overlap] += stamp[overlap]
        return subimage

class OverlapAnalyzer(object):
    """Analyze impact of overlapping sources on weak lensing.

    Args:
        survey(descwl.survey.Survey): Simulated survey to describe with FITS header keywords.
    """
    def __init__(self,survey):
        self.survey = survey
        self.models = [ ]
        self.stamps = [ ]
        self.bounds = [ ]

    def add_galaxy(self,model,stamps,bounds):
        """Add one galaxy to be analyzed.

        Args:
            model(:class:`descwl.model.Galaxy`): The galaxy model used for rendering.
            stamps(:class:`numpy.ndarray`): Array of shape (nstamp,width,height) containing
                pixel values for nstamp stamps of dimensions (width,height).
            bounds(galsim.BoundsI): Bounds of the stamps in the full simulated survey image.
        """
        self.models.append(model)
        self.stamps.append(stamps)
        self.bounds.append(bounds)

    def finalize(self):
        """Finalize analysis of all added galaxies.

        Returns:
            :class:`OverlapResults`: Overlap analysis results.
        """
        # Define columns and allocate space for our table data.
        num_galaxies = len(self.models)
        data = np.empty(num_galaxies,dtype=[
            ('db_id',np.int64),
            ('grp_id',np.int16),
            ('grp_size',np.int16),
            ('visible',np.bool8),
            ('dx',np.float32),
            ('dy',np.float32),
            ('z',np.float32),
            ('ab_mag',np.float32),
            ('snr_sky',np.float32),
            ('snr_iso',np.float32),
            ])

        # Calculate isolated galaxy quantities and identify overlapping groups.
        data['grp_id'] = np.arange(num_galaxies)
        for index,(model,stamps,bounds) in enumerate(zip(self.models,self.stamps,self.bounds)):
            data['db_id'][index] = model.identifier
            data['dx'][index] = model.dx_arcsecs
            data['dy'][index] = model.dy_arcsecs
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
            for pre_index,(pre_model,pre_stamps,pre_bounds) in enumerate(
                zip(self.models[:index],self.stamps[:index],self.bounds[:index])):
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

        table = astropy.table.Table(data,copy = False)
        results = OverlapResults(self.survey,table,self.stamps,self.bounds)
        return results
