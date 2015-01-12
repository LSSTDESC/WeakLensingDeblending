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
        if self.table is not None:
            self.num_objects = len(self.table)
            self.locals = { name: self.table[name] for name in self.table.colnames }
        self.stamps = stamps
        self.bounds = bounds
        if len(self.stamps) > 0:
            self.num_slices = self.stamps[0].shape[0]
            if self.num_slices not in (1,len(self.slice_labels)):
                raise RuntimeError('Image datacubes have unexpected number of slices (%d).'
                    % self.num_slices)

    slice_labels = ['dflux','dx','dy','dscale','dg1','dg2']

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
            galsim.ImageView: A GalSim image for the requested galaxy. Note that the
                return value is a read-only view of our internal stamps array.

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

    def get_fisher_images(self,index1,index2,background):
        """Return Fisher-matrix images for a pair of galaxies.

        Fisher matrix images are derived from the partial derivatives with respect to
        six parameters: total flux in detected electrons, centroid positions in x and y
        in arcseconds, flux-preserving radial scale (dimensionless), and shear g1,g2
        with \|g\| = (a-b)/(a+b). Use the :ref:`prog-fisher` program to display Fisher
        matrix images.

        Args:
            index1(int): Index of the first galaxy to use.
            index2(int): Index of the second galaxy to use, which might the same as index1.
            background(galsim.Image): Background image that combines all sources that overlap
                the two galaxies and completely contains them.

        Returns:
            tuple: Tuple (images,overlap) where images is a :type:`numpy.ndarray` with shape
                (npar,npar,height,width), where npar=6 and (height,width) are the dimensions
                of the overlap between the two galaxies, and overlap gives the bounding box of the
                overlap in the full survey image. Returns None,None if the two galaxies do not overlap.

        Raises:
            RuntimeError: Invalid index1 or index2, or galaxies are not contained with the
                background image, or no partial derivative images are available.
        """
        npar = self.num_slices
        if npar != len(self.slice_labels):
            raise RuntimeError('No partial derivative images are available.')
        # Calculate the overlap bounds.
        try:
            overlap = self.bounds[index1] & self.bounds[index2]
        except IndexError:
            raise RuntimeError('Invalid index1=%d or index2=%d.' % (index1,index2))
        # Check that the background image contains each galaxy.
        if not background.bounds.includes(self.bounds[index1]):
            raise RuntimeError('Galaxy %d is not contained within the background image.' % index1)
        if not background.bounds.includes(self.bounds[index2]):
            raise RuntimeError('Galaxy %d is not contained within the background image.' % index2)
        # Is there any overlap between the two galaxies?
        if overlap.area() == 0:
            return None,None
        product = self.get_stamp(index1,0)[overlap]*self.get_stamp(index2,0)[overlap]
        if not np.any(product.array):
            return None,None
        # Fill arrays of partial derivatives within the overlap region.
        width = overlap.xmax - overlap.xmin + 1
        height = overlap.ymax - overlap.ymin + 1
        partials1 = np.empty((npar,height,width),dtype = np.float32)
        partials2 = np.empty((npar,height,width),dtype = np.float32)
        for islice in range(npar):
            partials1[islice] = self.get_stamp(index1,islice)[overlap].array
            partials2[islice] = self.get_stamp(index2,islice)[overlap].array
        partials1[0] /= self.table['flux'][index1]
        partials2[0] /= self.table['flux'][index2]
        # Normalize the Fisher images.
        mu0 = background[overlap].array + self.survey.mean_sky_level
        fisher_norm = mu0**-1 + 0.5*mu0**-2
        # Calculate the Fisher images as the outer product of the partial-derivative images.
        images = np.einsum('yx,iyx,jyx->ijyx',fisher_norm,partials1,partials2)
        return images,overlap

    def get_matrices(self,selected):
        """Return matrices derived the from Fisher-matrix images for a set of sources.

        If the Fisher matrix is not invertible, `covariance`, `variance` and `correlation`
        will be returned as None.  If any variances are <= 0, `correlation` will be
        returned as None. These problems are generally associated with sources that are
        barely above the pixel SNR threshold. Use the :ref:`prog-fisher` program to
        visualize matrix elements and to analyze cases where some arrays are returned as None.

        Args:
            selected(iterable): Array of integer indices for the sources to include in the
                calculated matrices.

        Returns:
            tuple: Tuple `(fisher,covariance,variance,correlation)` of :class:`numpy.ndarray`
                where `variance` has shape (npar,) and all other arrays are symmetric with
                shape (npar,npar), where npar = 6*len(selected). If any array cannot be
                calculated, it will be returned as None.
        """
        background = self.get_subimage(selected)
        npar = self.num_slices
        nfisher = len(selected)*npar
        fisher = np.zeros((nfisher,nfisher),dtype = np.float64)
        for row,index1 in enumerate(selected):
            for col,index2 in enumerate(selected[:row+1]):
                images,overlap = self.get_fisher_images(index1,index2,background)
                if overlap is None:
                    continue
                fisher_sums = np.sum(images,axis=(2,3),dtype = np.float64)
                fisher[npar*row:npar*(row+1),npar*col:npar*(col+1)] = fisher_sums
                if row != col:
                    fisher[npar*col:npar*(col+1),npar*row:npar*(row+1)] = fisher_sums.T

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

    def finalize(self,verbose,trace):
        """Finalize analysis of all added galaxies.

        Args:
            verbose(bool): Print a summary of analysis results.
            trace(callable): Function to call for tracing resource usage. Will be
                called with a brief :type:`str` description of each checkpoint.

        Returns:
            :class:`OverlapResults`: Overlap analysis results.
        """
        trace('OverlapAnalyzer.finalize begin')
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
            ('snr_iso2',np.float32),
            ('snr_isof',np.float32),
            ('snr_grpf',np.float32),
            ])
        trace('allocated table of %ld bytes for %d galaxies' % (data.nbytes,num_galaxies))

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
        data['grp_id'] = np.arange(num_galaxies)
        for index in range(num_galaxies):
            trace('index %d' % index)
            model,stamps,bounds = self.models[index],self.stamps[index],self.bounds[index]
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
            print 'Simulated %d galaxies in %d overlap groups.' % (num_galaxies,num_groups)

        # Initialize our results object so we can use its methods (but be careful not
        # to use a method that needs something in table that we have not filled in yet).
        table = astropy.table.Table(data,copy = False)
        results = OverlapResults(self.survey,table,self.stamps,self.bounds)
        npartials = results.num_slices
        assert npartials == len(results.slice_labels), (
            'Datacubes have wrong num_slices = %d' % npartials)

        sky = self.survey.mean_sky_level
        # Loop over groups to calculate pixel-level quantities.
        for i,grp_id in enumerate(grp_id_set):
            trace('grp_id %d is %d of %d' % (grp_id,i,len(grp_id_set)))
            grp_members = (data['grp_id'] == grp_id)
            grp_size = np.count_nonzero(grp_members)
            data['grp_size'][grp_members] = grp_size
            group_indices = np.arange(num_galaxies)[grp_members]
            group_image = results.get_subimage(group_indices)
            fisher,covariance,variance,correlation = results.get_matrices(group_indices)
            for index,galaxy in enumerate(group_indices):
                flux = data['flux'][galaxy]
                signal = results.get_stamp(galaxy)
                # Calculate this galaxy's purity.
                signal_plus_background = group_image[signal.bounds]
                data['purity'][galaxy] = (
                    np.sum(signal.array**2)/np.sum(signal.array*signal_plus_background.array))
                # Calculate the SNR this galaxy would have without any overlaps and
                # assuming that we are in the sky-dominated limit.
                data['snr_sky'][galaxy] = np.sqrt(np.sum(signal.array**2)/sky)
                # Calculate this galaxy's SNR in various ways.
                flux_index = index*npartials
                data['snr_iso2'][galaxy] = flux*np.sqrt(fisher[flux_index,flux_index])
                if correlation is not None:
                    data['snr_grpf'][galaxy] = flux/np.sqrt(variance[flux_index])
                else:
                    data['snr_grpf'][galaxy] = -1.
                if grp_size == 1:
                    data['snr_iso'][galaxy] = data['snr_iso2'][galaxy]
                    data['snr_isof'][galaxy] = data['snr_grpf'][galaxy]
                else:
                    iso_fisher,iso_covariance,iso_variance,iso_correlation = (
                        results.get_matrices([galaxy]))
                    data['snr_iso'][galaxy] = flux*np.sqrt(iso_fisher[0,0])
                    if iso_correlation is not None:
                        data['snr_isof'][galaxy] = flux/np.sqrt(iso_variance[0])
                    else:
                        data['snr_isof'][galaxy] = -1.

            # Order group members by decreasing isolated S/N.
            sorted_indices = group_indices[np.argsort(data['snr_iso'][grp_members])[::-1]]
            data['grp_rank'][sorted_indices] = np.arange(grp_size,dtype = np.int16)
            # Replace group ID with ID of galaxy with largest S/N.
            group_leader = data['db_id'][sorted_indices[0]]
            data['grp_id'][grp_members] = group_leader

        trace('OverlapAnalyzer.finalize end')
        return results
