"""Perform weak-lensing analysis of simulated sources.
"""

import numpy as np

import astropy.table

import galsim

import descwl.model

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

    Raises:
        RuntimeError: Image datacubes have unexpected number of slices.
    """
    def __init__(self,survey,table,stamps,bounds):
        self.survey = survey
        self.table = table
        self.stamps = stamps
        self.bounds = bounds
        self.num_objects = len(self.table)
        self.num_slices = self.stamps[0].shape[0]
        self.slice_labels = ['dflux','dx','dy','dscale','dg1','dg2']
        if self.num_slices not in (1,len(self.slice_labels)):
            raise RuntimeError('Image datacubes have unexpected number of slices (%d).'
                % self.num_slices)
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

    def get_subimage(self,indices,datacube_index=0):
        """Return simulated subimage of a set of objects.

        Args:
            indices(iterable): Indices of the objects to include in the subimage.

        Returns:
            galsim.Image: Image of the selected objects or None if indices is empty.

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
            stamp = self.get_stamp(index,datacube_index)
            overlap = subimage_bounds & stamp.bounds
            subimage[overlap] += stamp[overlap]
        return subimage

    def get_fisher_images(self,indices):
        """Return Fisher-matrix images for a set of objects.

        Fisher matrix images are derived from the partial derivatives with respect to
        six parameters: total flux in detected electrons, centroid positions in x and y
        in arcseconds, flux-preserving radial scale (dimensionless), and shear g1,g2
        with \|g\| = (a-b)/(a+b).

        Args:
            indices(iterable): Indices of the objects to include in the subimage.

        Returns:
            numpy.ndarray: Array with shape (npar,npar,height,width) where npar = 6*len(indices)
                is the number of floating parameters used for the calculation, and (height,width)
                are the dimensions of the subimage containing the requested objects.
                The returned array is symmetric in its first two indices. Returns None if
                indices is empty.

        Raises:
            RuntimeError: An index is out of range.
        """
        nselected = len(indices)
        npartials = self.num_slices
        npar = nselected*npartials
        background = self.get_subimage(indices)
        if nselected == 0 or background is None:
            return None
        height,width = background.array.shape
        sky_level = self.survey.mean_sky_level

        fisher_images = np.empty((npar,npar,height,width))
        stamp1 = background.copy()
        stamp2 = background.copy()

        for index1 in range(npar):
            galaxy1 = indices[index1//npartials]
            slice1 = index1%npartials
            stamp1.array[:] = 0.
            stamp1[self.bounds[galaxy1]] = self.get_stamp(galaxy1,slice1)
            if slice1 == 0:
                # Normalize to give partial with respect to added flux in electrons.
                stamp1 /= self.table['flux'][galaxy1]
            for index2 in range(index1+1):
                galaxy2 = indices[index2//npartials]
                slice2 = index2%npartials
                stamp2.array[:] = 0.
                stamp2[self.bounds[galaxy2]] = self.get_stamp(galaxy2,slice2)
                if slice2 == 0:
                    # Normalize to give partial with respect to added flux in electrons.
                    stamp2 /= self.table['flux'][galaxy2]
                fisher_images[index1,index2] = stamp1.array*stamp2.array/(
                    background.array + sky_level)
                if index2 < index1:
                    fisher_images[index2,index1] = fisher_images[index1,index2]
        return fisher_images

    def get_matrices(self,fisher_images):
        """Return matrices derived from Fisher-matrix images.

        If the Fisher matrix is not invertible, `covariance`, `variance` and `correlation`
        will be returned as None.  If any variances are <= 0, `correlation` will be
        returned as None.

        Args:
            fisher_images(numpy.ndarray): Fisher-matrix images returned by
                :meth:`get_fisher_images`. These are assumed to be symmetric in the
                first two indices, but this is not checked.

        Returns:
            tuple: Tuple `(fisher,covariance,variance,correlation)` of :class:`numpy.ndarray`
                where `variance` has shape (npar,) and all other arrays are symmetric with
                shape (npar,npar).  If any array cannot be calculated, it will be returned
                as None.

        Raises:
            RuntimeError: `fisher_images` has unexpected shape.
        """
        shape = fisher_images.shape
        if len(shape) != 4 or shape[0] != shape[1]:
            raise RuntimeError('Fisher images have unexpected shape %r' % shape)

        fisher = np.sum(fisher_images,axis=(2,3))
        covariance,variance,correlation = None,None,None
        try:
            covariance = np.linalg.inv(fisher)
            variance = np.diag(covariance)
            if np.min(variance) > 0:
                correlation = covariance/np.sqrt(np.outer(variance,variance))
        except np.linalg.LinAlgError:
            pass

        return fisher,covariance,variance,correlation

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

    def finalize(self,verbose):
        """Finalize analysis of all added galaxies.

        Args:
            verbose(bool): Print a summary of analysis results.

        Returns:
            :class:`OverlapResults`: Overlap analysis results.
        """
        # Define columns and allocate space for our table data.
        num_galaxies = len(self.models)
        data = np.empty(num_galaxies,dtype=[
            ('db_id',np.int64),
            ('grp_id',np.int64),
            ('grp_size',np.int16),
            ('grp_rank',np.int16),
            ('visible',np.bool8),
            # Source properties.
            ('f_disk', np.float32),
            ('f_bulge', np.float32),
            ('dx',np.float32),
            ('dy',np.float32),
            ('z',np.float32),
            ('ab_mag',np.float32),
            ('flux',np.float32),
            ('sigma_m',np.float32),
            ('sigma_p',np.float32),
            ('e1',np.float32),
            ('e2',np.float32),
            ('a',np.float32),
            ('b',np.float32),
            ('beta',np.float32),
            # Pixel-level properties.
            ('purity',np.float32),
            ('snr_sky',np.float32),
            ('snr_iso',np.float32),
            ('snr_grp',np.float32),
            ('snr_isof',np.float32),
            ('snr_grpf',np.float32),
            ])

        # Initialize integer arrays of bounding box limits.
        xmin = np.empty(num_galaxies,np.int32)
        ymin = np.empty(num_galaxies,np.int32)
        xmax = np.empty(num_galaxies,np.int32)
        ymax = np.empty(num_galaxies,np.int32)
        for i in range(num_galaxies):
            xmin[i] = self.bounds[i].xmin
            ymin[i] = self.bounds[i].ymin
            xmax[i] = self.bounds[i].xmax
            ymax[i] = self.bounds[i].ymax

        # Find overlapping bounding boxes.
        x_overlap = np.logical_and(
            xmin[:,np.newaxis] <= xmax[np.newaxis,:],
            ymin[:,np.newaxis] <= ymax[np.newaxis,:])
        y_overlap = np.logical_and(
            xmax[:,np.newaxis] >= xmin[np.newaxis,:],
            ymax[:,np.newaxis] >= ymin[np.newaxis,:])
        overlapping_bounds = np.logical_and(x_overlap,y_overlap)

        # Calculate isolated galaxy quantities and identify overlapping groups.
        # At this stage, we use small integers for grp_id values, but these are later
        # reset to be the db_id of each group's leader.
        sky = self.survey.mean_sky_level
        data['grp_id'] = np.arange(num_galaxies)
        for index,(model,stamps,bounds) in enumerate(zip(self.models,self.stamps,self.bounds)):
            data['db_id'][index] = model.identifier
            data['dx'][index] = model.dx_arcsecs
            data['dy'][index] = model.dy_arcsecs
            data['z'][index] = model.redshift
            data['ab_mag'][index] = model.ab_magnitude
            data['flux'][index] = model.model.getFlux()
            # Is this galaxy's centroid visible in the survey image?
            data['visible'][index] = self.survey.image.bounds.includes(bounds.center())
            # Save model parameters.
            data['f_disk'][index] = model.disk_fraction
            data['f_bulge'][index] = model.bulge_fraction
            # Calculate this galaxy's size and shape from its second-moments tensor.
            sigma_m,sigma_p,a,b,beta,e1,e2 = descwl.model.moments_size_and_shape(
                model.second_moments)
            # Save values to the analysis results.
            data['sigma_m'][index] = sigma_m
            data['sigma_p'][index] = sigma_p
            data['a'][index] = a
            data['b'][index] = b
            data['e1'][index] = e1
            data['e2'][index] = e2
            data['beta'][index] = beta
            # Calculate the SNR this galaxy would have without any overlaps in the
            # sky-dominated limit.
            fiducial = stamps[0].flatten()
            data['snr_sky'][index] = np.sqrt(np.sum(fiducial**2)/sky)
            # Include the signal in the noise variance.
            mu0 = fiducial + sky
            data['snr_iso'][index] = np.sqrt(np.sum(fiducial**2*(mu0**-1 + 0.5*mu0**-2)))
            # Loop over earlier galaxies with overlapping bounding boxes.
            for pre_index in np.arange(index)[overlapping_bounds[index,:index]]:
                pre_bounds = self.bounds[pre_index]
                overlap = bounds & pre_bounds
                assert overlap.area() > 0
                # Are there any overlapping pixels with non-zero flux from both sources?
                pre_stamps = self.stamps[pre_index]
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

        # How many groups do we have?
        grp_id_set = set(data['grp_id'])
        num_groups = len(grp_id_set)
        if verbose:
            print 'Simulated %d galaxies in %d overlap groups' % (num_galaxies,num_groups)

        # Initialize our results object so we can use its methods (but be careful not
        # to use a method that needs something in table that we have not filled in yet).
        table = astropy.table.Table(data,copy = False)
        results = OverlapResults(self.survey,table,self.stamps,self.bounds)
        npartials = results.num_slices
        assert npartials == len(results.slice_labels), (
            'Datacubes have wrong num_slices = %d' % npartials)

        # Initialize assuming that sources do not overlap.
        data['snr_grp'] = data['snr_iso']
        data['purity'] = 1.
        data['grp_rank'] = 0
        # Loop over groups to calculate quantities that need the proto-results.
        for grp_id in grp_id_set:
            grp_members = (data['grp_id'] == grp_id)
            grp_size = np.count_nonzero(grp_members)
            data['grp_size'][grp_members] = grp_size
            if grp_size > 1:
                group_indices = np.arange(num_galaxies)[grp_members]
                group_image = results.get_subimage(group_indices)
                # Loop over pairs of galaxies in this overlapping group to calculate
                # the Fisher matrix for the overlapping S/N calculation.
                fisher = np.zeros((grp_size,grp_size),dtype = np.float64)
                flux = np.empty(grp_size,dtype = np.float64)
                for i1,g1 in enumerate(group_indices):
                    stamp1 = results.get_stamp(g1)
                    flux[i1] = self.models[g1].model.getFlux()
                    # Calculate this source's purity within its group.
                    stamp1_group = group_image[stamp1.bounds]
                    data['purity'][g1] = (
                        np.sum(stamp1.array**2)/np.sum(stamp1.array*stamp1_group.array))
                    for i2,g2 in enumerate(group_indices[:i1+1]):
                        # Galaxies (g1,g2) might not directly overlap even if they are in
                        # an overlapping group.
                        if not overlapping_bounds[g1,g2]:
                            continue
                        stamp2 = results.get_stamp(g2)
                        overlap = stamp1.bounds & stamp2.bounds
                        assert overlap.area() > 0
                        fiducial1 = stamp1[overlap].array.flatten()
                        fiducial2 = stamp2[overlap].array.flatten()
                        mu0 = group_image[overlap].array.flatten() + sky
                        fisher[i1,i2] = np.sum(
                            fiducial1*fiducial2*(mu0**-1 + 0.5*mu0**-2))/(flux[i1]*flux[i2])
                        fisher[i2,i1] = fisher[i1,i2]
                covariance = np.linalg.inv(fisher)
                variance = np.diag(covariance)
                if not np.all(variance > 0):
                    print 'variance < 0',variance
                    print fisher
                group_snr = np.sqrt(flux**2/variance)
                data['snr_grp'][grp_members] = group_snr
                # Order group members by decreasing group S/N.
                sorted_indices = group_indices[np.argsort(group_snr)[::-1]]
                data['grp_rank'][sorted_indices] = np.arange(grp_size,dtype = np.int16)
                # Replace group ID with ID of galaxy with largest S/N.
                group_leader = data['db_id'][sorted_indices[0]]
                data['grp_id'][grp_members] = group_leader

        return results
