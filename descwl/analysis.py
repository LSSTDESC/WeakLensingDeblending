"""Perform weak-lensing analysis of simulated sources.
"""

import numpy as np
import scipy.spatial

import astropy.table

import galsim

import lmfit

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
        stamps(array): Array of :class:`numpy.ndarray` postage-stamp datacubes. Might be None.
        bounds(array): Array of galsim.BoundsI objects giving pixel coordinates of
            each datacube in the full simulated image. Might be None.
        num_slices(int): Number of datacube slices for each stamp. Ignored if stamps is None.

    Raises:
        RuntimeError: Image datacubes have unexpected number of slices.
    """
    def __init__(self,survey,table,stamps,bounds,num_slices):
        self.survey = survey
        self.table = table
        if self.table is not None:
            self.num_objects = len(self.table)
            self.locals = { name: self.table[name] for name in self.table.colnames }
        self.stamps = stamps
        self.bounds = bounds
        self.num_slices = num_slices
        if len(self.stamps) > 0:
            if self.num_slices not in (1,len(self.slice_labels)):
                raise RuntimeError('Image datacubes have unexpected number of slices (%d).'
                    % self.num_slices)
        self.noise_seed = None

    slice_labels = ['dflux','dx','dy','ds','dg1','dg2']

    def add_noise(self,noise_seed):
        """Add Poisson noise to the simulated survey image.

        Args:
            noise_seed(int): Random seed to use.

        Raises:
            RuntimeError: Noise has already been added.
        """
        if self.noise_seed is not None:
            raise RuntimeError('Noise has already been added to the simulated image.')
        self.noise_seed = noise_seed
        generator = galsim.random.BaseDeviate(seed = self.noise_seed)
        noise = galsim.PoissonNoise(rng = generator, sky_level = self.survey.mean_sky_level)
        self.survey.image.addNoise(noise)

    def select(self,*selectors,**options):
        """Select objects.

        This function is implemented using :func:`eval` so should only be used when the
        source of the selector strings is trusted.

        Args:
            selectors(str,...): A sequence of one or more selection strings, each
                consisting of a boolean expression of column table names, e.g.
                'snr_iso > 5'.  The special strings 'ALL' and 'NONE' do what you
                would expect.
            mode(str): The string 'and' or 'or' to specify how multiple selectors
                should be combined.
            format(str): The string 'mask' or 'index' to specify the format of the
                returned array.  A mask is an array of booleans with one entry per
                source. The alternative index array of integers lists the selected
                indices in increasing order and might be empty.  The mask format
                is useful for performing your own logical operations on the
                returned array. The index format is useful for iterating over the
                selected objects. Both formats can be used to slice the catalog
                table to extract only the selected rows.

        Returns:
            :class:`numpy.ndarray`: A numpy array of booleans or integers, depending
                on the value of the format input argument.

        Raises:
            RuntimeError: A selector uses an undefined variable name, the method was
                passed an unexpected mode or format string, or these results have
                no catalog table.
        """
        if self.table is None:
            raise RuntimeError('Cannot select() without a catalog table.')
        mode = options.get('mode','and')
        if mode not in ('and','or'):
            raise RuntimeError('Unexpected mode "%s".' % mode)
        format = options.get('format','index')
        if format not in ('mask','index'):
            raise RuntimeError('Unexpected format "%s".' % format)

        selection = self._select('ALL') if mode == 'and' else self._select('NONE')
        for selector in selectors:
            update = self._select(selector)
            if mode == 'and':
                selection = np.logical_and(selection,update)
            else:
                selection = np.logical_or(selection,update)

        if format == 'index':
            selection = np.arange(self.num_objects)[selection]
        return selection

    def _select(self,selector):
        # This private method implements select() for a single selector.
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
        datacube = self.stamps[index]
        if callable(datacube):
            datacube = datacube()
            self.stamps[index] = datacube
        image = galsim.Image(datacube[datacube_index],
            xmin = bounds.xmin,ymin = bounds.ymin,
            scale = self.survey.pixel_scale,make_const = False)
        return image

    def get_subimage(self,indices,datacube_index=0):
        """Return simulated subimage of a set of objects.

        Args:
            indices(iterable): Indices of the objects to include in the subimage.
                Our :meth:`select` method returns a suitable argument to use here
                by default (format = 'index').

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
            tuple: Tuple (images,overlap) where images is a :class:`numpy.ndarray` with shape
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

        If the Fisher matrix is not invertible or any variances are <= 0, we will drop
        the selected source with the lowest value of snr_iso and try again. This procedure
        is iterated until we get a valid covariance matrix. Matrix elements for all parameters
        of any sources that get dropped by this procedure will be set to zero and variances
        will be set to np.inf so that 1/np.sqrt(variance) = 0. Invalid covariances are
        generally associated with sources that are barely above the pixel SNR threshold,
        so this procedure should normally provide sensible values for the largest
        possible subset of the input selected sources. Use the :ref:`prog-fisher` program to
        visualize matrix elements and to further study examples of invalid covariances.

        Args:
            selected(iterable): Array of integer indices for the sources to include in the
                calculated matrices.

        Returns:
            tuple: Tuple `(fisher,covariance,variance,correlation)` of :class:`numpy.ndarray`
                where `variance` has shape (npar,) and all other arrays are symmetric with
                shape (npar,npar), where npar = 6*len(selected). Matrix elements will be
                zero for any parameters associated with dropped sources, as described above.
        """
        background = self.get_subimage(selected)
        nsel = len(selected)
        npar = self.num_slices
        nfisher = nsel*npar
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

        # Sort indices into the selected array by increasing snr_iso.
        priority = np.arange(nsel)[np.argsort(self.table['snr_iso'][selected])]
        # Start by trying to use all sources in the group.
        keep = np.ones((nsel,npar),dtype=bool)
        num_dropped = 0
        while np.any(keep):
            try:
                keep_flat = keep.flatten()
                # Advanced indexing like this makes a copy, not a view.
                reduced_fisher = fisher[keep_flat,:][:,keep_flat]
                reduced_covariance = np.linalg.inv(reduced_fisher)
                reduced_variance = np.diag(reduced_covariance)
                assert np.min(reduced_variance) > 0,'Expected variance > 0'
                reduced_correlation = reduced_covariance/np.sqrt(
                    np.outer(reduced_variance,reduced_variance))
                break
            except (np.linalg.LinAlgError,AssertionError),e:
                # We can't calculate a covariance for this set of objects, so drop the next
                # lowest SNR member of the set and try again.
                keep[priority[num_dropped],:] = False
                num_dropped += 1

        if num_dropped == 0:
            return reduced_fisher, reduced_covariance,reduced_variance,reduced_correlation
        else:
            fisher = np.zeros((nfisher,nfisher),dtype = np.float64)
            covariance = np.zeros((nfisher,nfisher),dtype = np.float64)
            correlation = np.zeros((nfisher,nfisher),dtype = np.float64)
            variance = np.empty((nfisher,),dtype = np.float64)
            variance[:] = np.inf
            # Build matrices with zeros for any sources that had to be dropped. Is there
            # a more elegant way to do this? We cannot simply assign to a submatrix
            # since that requires advancing slicing, which creates a copy, not a view.
            keep_map = np.zeros(nsel,dtype=int)
            keep_map[keep[:,0]] = np.arange(nsel-num_dropped)
            next = 0
            for row in range(nsel):
                if not keep[row,0]: continue
                row_slice = slice(npar*row,npar*(row+1))
                krow = keep_map[row]
                krow_slice = slice(npar*krow,npar*(krow+1))
                variance[row_slice] = reduced_variance[krow_slice]
                for col in range(row+1):
                    if not keep[col,0]: continue
                    col_slice = slice(npar*col,npar*(col+1))
                    kcol = keep_map[col]
                    kcol_slice = slice(npar*kcol,npar*(kcol+1))
                    fisher[row_slice,col_slice] = reduced_fisher[krow_slice,kcol_slice]
                    covariance[row_slice,col_slice] = reduced_covariance[krow_slice,kcol_slice]
                    correlation[row_slice,col_slice] = reduced_correlation[krow_slice,kcol_slice]
                    if row == col: continue
                    fisher[col_slice,row_slice] = reduced_fisher[kcol_slice,krow_slice]
                    covariance[col_slice,row_slice] = reduced_covariance[kcol_slice,krow_slice]
                    correlation[col_slice,row_slice] = reduced_correlation[kcol_slice,krow_slice]
            return fisher,covariance,variance,correlation

    def match_sextractor(self,catalog_name,column_name = 'match'):
        """Match detected objects to simulated sources.

        Args:
            catalog_name(str): Name of an ASCII catalog in SExtractor-compatible format, and containing
                X_IMAGE, Y_IMAGE columns and `num_found` rows.
            column_name(str): Name of a column in our `table` data member to fill with the detected
                object index matched to each simulated source, or -1 if no match is found. Overwrites
                any exisiting column or creates a new one.  Does nothing if column_name is None or ''.

        Returns:
            tuple: Tuple `detected,matched,indices,distance` where `detected` is the detected catalog as
                a :class:`astropy.table.Table`, `matched` is an array of `num_found` booleans
                indicating whether each object was matched with a simulated object, `indices` is an array
                of `num_matched = np.count_nonzero(matched)` integers giving row numbers in our `table`
                attribute for each simulated match, and `distance` is an array of `num_matched` separation
                distances in arcseconds between matched and simulated objects.
        """
        # Read the catalog of detected objects.
        detected = astropy.table.Table.read(catalog_name,format='ascii')
        assert 'X_IMAGE' in detected.colnames, 'Missing required column X_IMAGE'
        assert 'Y_IMAGE' in detected.colnames, 'Missing required column Y_IMAGE'
        # Build a KD tree for the simulated source positions (in arcsecs) relative to the image center.
        num_truth = len(self.table)
        xy_truth = np.empty((num_truth,2))
        xy_truth[:,0] = self.table['dx']
        xy_truth[:,1] = self.table['dy']
        kdtree = scipy.spatial.cKDTree(xy_truth)
        # Convert detected object centroids from pixels relative to bottom-left corner to arcsecs
        # into arcsecs relative to the image center.
        num_found = len(detected)
        xy_found = np.empty((num_found,2))
        width,height = self.survey.image_width,self.survey.image_height
        xy_found[:,0] = (detected['X_IMAGE'] - 0.5*width - 0.5)*self.survey.pixel_scale
        xy_found[:,1] = (detected['Y_IMAGE'] - 0.5*height - 0.5)*self.survey.pixel_scale
        # Find nearest simulated source to each detected source.
        min_distance,truth_indices = kdtree.query(xy_found)
        matched = (truth_indices < num_truth)
        indices = truth_indices[matched]
        distance = min_distance[matched]        
        # Add a table column with indices of matched detected objects or -1 for unmatched sources.
        if column_name:
            match_lookup = np.empty(len(self.table),dtype = int)
            match_lookup[:] = -1
            match_lookup[indices] = np.arange(len(detected))[matched]
            self.table[column_name] = match_lookup
            # Make this column available via our locals dictionary.
            self.locals[column_name] = self.table[column_name]
        return detected,matched,indices,distance

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

    def fit_galaxies(self,indices,observed_image,fixed_parameters = None):
        """Simultaneously fit a set of galaxy parameters to an observed image.

        Fits are performed on noise-free images, so there are no meaningful errors to report
        and only best-fit parameter values are returned.  This method is intended to support
        systematics studies where shifts in best-fit parameters are the quantities of interest.
        See the :meth:`descwl.render.GalaxyRenderer.draw` method for details on how the fit
        parameters are defined and how each galaxy model is rendered as these parameters
        are varied.

        Args:
            indices(iterable): List of num_galaxies integer galaxy indices to include in the fit.
                Index values correspond to the order in which galaxies were added using
                :meth:`add_galaxy`.
            observed_image(galsim.Image): A GalSim image that defines the observed pixels
                that will be fit and the region into which galaxy models will be renderered.
            fixed_parameters(dict): An optional dictionary of parameter values to fix in the
                fit, or None if all parameters should be floating.

        Returns:
            numpy.ndarray: Array with shape (num_galaxies,6) containing the best-fit values
                of the following six parameters: df,dx,dy,ds,dg1,dg2. See the
                :meth:`descwl.render.GalaxyRenderer.draw` method for details on how these
                parameters are defined.

        Raises:
            RuntimeError: Invalid galaxy index or fit did not find a minimum.
        """
        num_galaxies = len(indices)
        # Define the fit parameters we will use.
        parameters = lmfit.Parameters()
        for i in range(num_galaxies):
            parameters.add('df_%d' % i,value = 0.)
            parameters.add('dx_%d' % i,value = 0.)
            parameters.add('dy_%d' % i,value = 0.)
            parameters.add('ds_%d' % i,value = 0.)
            parameters.add('dg1_%d' % i,value = 0.,min = -0.2,max = +0.2)
            parameters.add('dg2_%d' % i,value = 0.,min = -0.2,max = +0.2)
        # Fix parameters as requested.
        if fixed_parameters:
            for name,value in fixed_parameters.iteritems():
                parameters[name].value = value
                parameters[name].vary = False
        # Initialize rendering.
        model_image = observed_image.copy()
        overlap = [ ]
        for i,galaxy in enumerate(indices):
            try:
                bbox = self.bounds[galaxy]
            except (TypeError,ValueError):
                raise RuntimeError('Invalid galaxy index in fit_galaxies: %r' % galaxy)
            overlap.append(model_image.bounds & bbox)
        data = observed_image.array.flatten()
        sigmas = np.sqrt(data + self.survey.mean_sky_level)
        # Define a function that evaluates the pixel residuals to be minimized.
        def residuals(p):
            model_image.array[:] = 0.
            for i,galaxy in enumerate(indices):
                model_stamp = self.models[galaxy].renderer.draw(df = p['df_%d'%i].value,
                    dx = p['dx_%d'%i].value,dy = p['dy_%d'%i].value,ds = p['ds_%d'%i].value,
                    dg1 = p['dg1_%d'%i].value,dg2 = p['dg2_%d'%i].value)
                model_image[overlap[i]] += model_stamp[overlap[i]]
            return (data - model_image.array.flat)/sigmas
        # Do the minimization.
        minimum = lmfit.minimize(residuals,parameters)
        if not minimum.success:
            raise RuntimeError('fit_galaxies did not find a minimum.')
        # Copy the best-fit parameter values into the returned array.
        bestfit_values = np.empty((num_galaxies,6))
        for i in range(num_galaxies):
            bestfit_values[i,0] = parameters['df_%d'%i].value
            bestfit_values[i,1] = parameters['dx_%d'%i].value
            bestfit_values[i,2] = parameters['dy_%d'%i].value
            bestfit_values[i,3] = parameters['ds_%d'%i].value
            bestfit_values[i,4] = parameters['dg1_%d'%i].value
            bestfit_values[i,5] = parameters['dg2_%d'%i].value
        return bestfit_values

    def finalize(self,verbose,trace):
        """Finalize analysis of all added galaxies.

        Args:
            verbose(bool): Print a summary of analysis results.
            trace(callable): Function to call for tracing resource usage. Will be
                called with a brief :class:`str` description of each checkpoint.

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
            ('visible',np.int16),
            # Stamp bounding box.
            ('xmin',np.int32),
            ('xmax',np.int32),
            ('ymin',np.int32),
            ('ymax',np.int32),
            # Source properties.
            ('f_disk', np.float32),
            ('f_bulge', np.float32),
            ('dx',np.float32),
            ('dy',np.float32),
            ('z',np.float32),
            ('ab_mag',np.float32),
            ('ri_color',np.float32),
            ('flux',np.float32),
            ('sigma_m',np.float32),
            ('sigma_p',np.float32),
            ('e1',np.float32),
            ('e2',np.float32),
            ('a',np.float32),
            ('b',np.float32),
            ('beta',np.float32),
            ('psf_sigm',np.float32),
            # Pixel-level properties.
            ('purity',np.float32),
            ('snr_sky',np.float32),
            ('snr_iso',np.float32),
            ('snr_grp',np.float32),
            ('snr_isof',np.float32),
            ('snr_grpf',np.float32),
            ('ds',np.float32),
            ('dg1',np.float32),
            ('dg2',np.float32),
            ('ds_grp',np.float32),
            ('dg1_grp',np.float32),
            ('dg2_grp',np.float32),
            # HSM analysis results.
            ('hsm_sigm',np.float32),
            ('hsm_e1',np.float32),
            ('hsm_e2',np.float32),
            # Systematics fit results.
            ('g1_fit',np.float32),
            ('g2_fit',np.float32),
            ])
        trace('allocated table of %ld bytes for %d galaxies' % (data.nbytes,num_galaxies))

        # Initialize integer arrays of bounding box limits.
        for i in range(num_galaxies):
            data['xmin'][i] = self.bounds[i].xmin
            data['ymin'][i] = self.bounds[i].ymin
            data['xmax'][i] = self.bounds[i].xmax
            data['ymax'][i] = self.bounds[i].ymax

        # Find overlapping bounding boxes.
        x_overlap = np.logical_and(
            data['xmin'][:,np.newaxis] <= data['xmax'][np.newaxis,:],
            data['ymin'][:,np.newaxis] <= data['ymax'][np.newaxis,:])
        y_overlap = np.logical_and(
            data['xmax'][:,np.newaxis] >= data['xmin'][np.newaxis,:],
            data['ymax'][:,np.newaxis] >= data['ymin'][np.newaxis,:])
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
            data['ri_color'][index] = model.ri_color
            data['flux'][index] = model.model.getFlux()
            # Is this galaxy's centroid visible in the survey image?
            data['visible'][index] = 1 if self.survey.image.bounds.includes(bounds.center()) else 0
            # Save model parameters.
            data['f_disk'][index] = model.disk_fraction
            data['f_bulge'][index] = model.bulge_fraction
            # Calculate this galaxy's sizes and shapes from its second-moments tensor.
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
            # Re-calculate sizes and shapes with the PSF second moments added.
            sigma_m_psf,sigma_p_psf,a_psf,b_psf,beta_psf,e1_psf,e2_psf = descwl.model.moments_size_and_shape(
                model.second_moments + self.survey.psf_second_moments)
            # Save the PSF-convolved sigma(-) since this can be directly compared with the HSM size.
            data['psf_sigm'][index] = sigma_m_psf
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
        num_slices,h,w = self.stamps[0].shape
        results = OverlapResults(self.survey,table,self.stamps,self.bounds,num_slices)

        # Check that we have partial derivatives available.
        if num_slices != len(results.slice_labels):
            raise RuntimeError('Missing required partial derivative images for Fisher matrix analysis.')

        sky = self.survey.mean_sky_level
        dflux_index = results.slice_labels.index('dflux')
        ds_index = results.slice_labels.index('ds')
        dg1_index = results.slice_labels.index('dg1')
        dg2_index = results.slice_labels.index('dg2')
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
                # Run the HSM analysis on this galaxy's stamp (ignoring overlaps).
                data['hsm_sigm'][galaxy] = np.nan
                data['hsm_e1'][galaxy] = np.nan
                data['hsm_e2'][galaxy] = np.nan
                try:
                    hsm_results = galsim.hsm.EstimateShear(signal,self.survey.psf_image)
                    data['hsm_sigm'][galaxy] = hsm_results.moments_sigma*self.survey.pixel_scale
                    data['hsm_e1'][galaxy] = hsm_results.corrected_e1
                    data['hsm_e2'][galaxy] = hsm_results.corrected_e2
                except RuntimeError,e:
                    # Usually "Unphysical situation: galaxy convolved with PSF is smaller than PSF!"
                    # due to truncation of a faint galaxy at the limiting isophote.  Try to just
                    # calculate the PSF-convolved size in this case.
                    try:
                        hsm_results = galsim.hsm.FindAdaptiveMom(signal)
                        data['hsm_sigm'][galaxy] = hsm_results.moments_sigma*self.survey.pixel_scale
                    except RuntimeError,e:
                        print str(e)
                # Calculate the SNR this galaxy would have without any overlaps and
                # assuming that we are in the sky-dominated limit.
                data['snr_sky'][galaxy] = np.sqrt(np.sum(signal.array**2)/sky)
                # Calculate this galaxy's SNR in various ways.
                base = index*num_slices
                data['snr_grp'][galaxy] = flux*np.sqrt(fisher[base+dflux_index,base+dflux_index])
                # Variances will be np.inf if this galaxy was dropped from the group for the
                # covariance calculation, leading to snr_grpf = 0 and infinite errors on s,g1,g2.
                data['snr_grpf'][galaxy] = flux/np.sqrt(variance[base+dflux_index])
                data['ds_grp'][galaxy] = np.sqrt(variance[base+ds_index])
                data['dg1_grp'][galaxy] = np.sqrt(variance[base+dg1_index])
                data['dg2_grp'][galaxy] = np.sqrt(variance[base+dg2_index])
                if grp_size == 1:
                    data['snr_iso'][galaxy] = data['snr_grp'][galaxy]
                    data['snr_isof'][galaxy] = data['snr_grpf'][galaxy]
                    data['ds'][galaxy] = data['ds_grp'][galaxy]
                    data['dg1'][galaxy] = data['dg1_grp'][galaxy]
                    data['dg2'][galaxy] = data['dg2_grp'][galaxy]
                else:
                    # Redo the Fisher matrix analysis but ignoring overlapping sources.
                    iso_fisher,iso_covariance,iso_variance,iso_correlation = (
                        results.get_matrices([galaxy]))
                    # snr_iso and snr_isof will be zero if the Fisher matrix is not invertible or
                    # yields any negative variances. Errors on s,g1,g2 will be np.inf.
                    data['snr_iso'][galaxy] = flux*np.sqrt(iso_fisher[dflux_index,dflux_index])
                    data['snr_isof'][galaxy] = flux/np.sqrt(iso_variance[dflux_index])
                    data['ds'][galaxy] = np.sqrt(iso_variance[ds_index])
                    data['dg1'][galaxy] = np.sqrt(iso_variance[dg1_index])
                    data['dg2'][galaxy] = np.sqrt(iso_variance[dg2_index])

            # Order group members by decreasing isolated S/N.
            sorted_indices = group_indices[np.argsort(data['snr_iso'][grp_members])[::-1]]
            data['grp_rank'][sorted_indices] = np.arange(grp_size,dtype = np.int16)
            # Replace group ID with ID of galaxy with largest S/N.
            group_leader = data['db_id'][sorted_indices[0]]
            data['grp_id'][grp_members] = group_leader

            alpha = +1.0
            data['g1_fit'][sorted_indices] = 0.
            data['g2_fit'][sorted_indices] = 0.
            if grp_size > 1:
                print '-- group size',grp_size
                # Loop over galaxies in order of decreasing snr_iso.
                for i1,g1 in enumerate(sorted_indices):
                    stamp1 = results.get_stamp(g1)
                    deblended = stamp1.copy()
                    # Loop over other galaxies in this group in order of decreasing snr_iso.
                    for i2,g2 in enumerate(sorted_indices):
                        if i1 == i2:
                            continue
                        stamp2 = results.get_stamp(g2)
                        bbox = stamp1.bounds & stamp2.bounds
                        if bbox.area() == 0:
                            continue
                        overlap = stamp1[bbox]*stamp2[bbox]/group_image[bbox]
                        # Re-assign a fraction of the overlapping flux in the deblended image.
                        if i1 < i2:
                            deblended[bbox] -= alpha*overlap
                        else:
                            deblended[bbox] += alpha*overlap
                    # Fit the deblended image of this galaxy.
                    try:
                        bestfit = self.fit_galaxies([g1],deblended)
                        print g1,data['db_id'][g1],bestfit
                        data['g1_fit'][g1] = bestfit[0,4]
                        data['g2_fit'][g1] = bestfit[0,5]
                    except RuntimeError,e:
                        print str(e)

        trace('OverlapAnalyzer.finalize end')
        return results
