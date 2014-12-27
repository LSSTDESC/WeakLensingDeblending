"""Render source models as simulated survey observations.
"""

import math
import inspect

import numpy as np

import galsim

class SourceNotVisible(Exception):
    """Custom exception to indicate that a source has no visible pixels above threshold.
    """
    pass

class Engine(object):
    """Rendering engine to simulate survey observations.

    Any pixels outside of the truncation radius or below the minimum S/N cut will have their
    flux set to zero in the rendered image. As a result the total rendered flux may be below
    the total model flux.

    Args:
        survey(descwl.survey.Survey): Survey that rendered images will simulate.
        min_snr(float): Simulate signals from individual sources down to this S/N threshold,
            where the signal N is calculated for the full exposure time and the noise N is
            set by the expected fluctuations in the sky background during a full exposure.
        truncate_radius(float): All extended sources are truncated at this radius in arcseconds.
        no_margin(bool): Do not simulate the tails of objects just outside the field.
        verbose_render(bool): Provide verbose output on rendering process.
    """
    def __init__(self,survey,min_snr,truncate_radius,no_margin,verbose_render):
        self.survey = survey
        self.min_snr = min_snr
        self.truncate_radius = truncate_radius
        self.no_margin = no_margin
        self.verbose_render = verbose_render
        # Calculate pixel flux threshold in electrons per pixel that determines how big a
        # bounding box we simulate for each source.
        sky_noise = math.sqrt(survey.mean_sky_level)
        self.pixel_cut = self.min_snr*sky_noise
        # Initialize our GalSim parameters.
        self.galsim_params = galsim.GSParams(maximum_fft_size=32768)
        # Evaluate the PSF dilution factor as the maximum fraction of a source's total flux
        # that can end up in a single pixel after convolution with the PSF.
        psf_stamp = galsim.ImageD(1,1,scale=self.survey.pixel_scale)
        self.survey.psf_model.drawImage(image = psf_stamp)
        self.psf_dilution = psf_stamp.array[0]
        # We will render each source into a square stamp with width = height = 2*padding + 1.
        self.padding = int(math.ceil(self.truncate_radius/self.survey.pixel_scale - 0.5))
        # Prepare a truncation mask.
        pixel_grid = np.arange(-self.padding,self.padding+1)*self.survey.pixel_scale
        pixel_x,pixel_y = np.meshgrid(pixel_grid,pixel_grid)
        pixel_radius = np.sqrt(pixel_x**2 + pixel_y**2)
        self.truncation_mask = (pixel_radius <= self.truncate_radius)

    def description(self):
        """Describe our rendering configuration.

        Returns:
            str: Description of the rendering configuration that we be used to simulate
                the survey.
        """
        return '\n'.join([
            ('Will render all pixels with at least %.1f detected electrons.' % self.pixel_cut),
            ('PSF dilution factor is %.6f.' % self.psf_dilution)
            ])

    def render_galaxy(self,galaxy):
        """Render a galaxy model for a simulated survey.

        Args:
            galaxy(descwl.model.Galaxy): Model of the galaxy to render.

        Returns:
            tuple: `(stamps,bounds)` where `stamps` is a :class:`numpy.ndarray` of shape
                (nstamp,width,height) pixel values that represents nstamp postage-stamp images
                with the same dimensions (width,height) determined by the rendering options
                provided. The returned `bounds` give the position of these stamps within the
                full simulated survey image as a `galsim.BoundsI` object. Note that these bounds
                might extend beyond the survey image, but will always have some overlap where
                the source is above threshold.

        Raises:
            SourceNotVisible: Galaxy has no pixels above threshold that are visible in the
                simulated survey.
        """
        # Skip sources that are too faint to possibly be above our cut after PSF convolution.
        if galaxy.model.getFlux()*self.psf_dilution < self.pixel_cut:
            raise SourceNotVisible

        # Calculate the offset of the source center from the bottom-left corner of the
        # simulated image in floating-point pixel units.
        centroid = galaxy.model.centroid()
        x_center_pixels,y_center_pixels = self.survey.get_image_coordinates(centroid.x,centroid.y)

        # Calculate the corresponding central pixel indices in the full image, where (0,0) is the
        # bottom-left corner.
        x_center_index = int(math.floor(x_center_pixels))
        y_center_index = int(math.floor(y_center_pixels))

        # Calculate the bounding box to use for simulating this galaxy.
        x_min = x_center_index - self.padding
        x_max = x_center_index + self.padding
        y_min = y_center_index - self.padding
        y_max = y_center_index + self.padding
        bounds = galsim.BoundsI(x_min,x_max,y_min,y_max)

        # Calculate the offset of the bounding box center from the image center in arcsecs.
        dx_stamp_arcsec = 0.5*(x_min + x_max+1 - self.survey.image_width)*self.survey.pixel_scale
        dy_stamp_arcsec = 0.5*(y_min + y_max+1 - self.survey.image_height)*self.survey.pixel_scale

        # Shift the model to the bounding box center and convolve with the survey PSF.
        # We do not convolve by the pixel response since drawImage takes care of this.
        model = galsim.Convolve([
            galaxy.model.shift(dx=-dx_stamp_arcsec,dy=-dy_stamp_arcsec),
            self.survey.psf_model
            ],gsparams=self.galsim_params)

        # Render the model in its own postage stamp.
        stamp = galsim.Image(bounds = bounds,scale = self.survey.pixel_scale, dtype = np.float32)
        model.drawImage(image = stamp, use_true_center = True)

        # Identify pixels with flux above our cut and within our truncation radius
        # and zero all other pixel fluxes.
        keep_mask = (stamp.array*self.truncation_mask > self.pixel_cut)
        if np.sum(keep_mask) == 0:
            raise SourceNotVisible
        stamp.array[np.logical_not(keep_mask)] = 0.

        # Crop the bounding box.
        x_projection = (np.sum(keep_mask,axis=0) > 0)
        y_projection = (np.sum(keep_mask,axis=1) > 0)
        x_min_inset = np.argmax(x_projection)
        x_max_inset = np.argmax(x_projection[::-1])
        y_min_inset = np.argmax(y_projection)
        y_max_inset = np.argmax(y_projection[::-1])
        cropped_bounds = galsim.BoundsI(
            x_min+x_min_inset,x_max-x_max_inset,
            y_min+y_min_inset,y_max-y_max_inset)
        cropped_stamp = stamp[cropped_bounds]

        # Add the rendered model to the survey image.
        survey_overlap = cropped_bounds & self.survey.image.bounds
        if survey_overlap.area() == 0:
            raise SourceNotVisible
        self.survey.image[survey_overlap] += cropped_stamp[survey_overlap]

        # Draw directly into the survey image, with no truncation.
        #model.drawImage(image = self.survey.image,add_to_image = True,use_true_center = True)

        if self.verbose_render:
            print 'Rendered galaxy model for id = %d with z = %.3f' % (
                galaxy.identifier,galaxy.redshift)
            print 'bounds: [%d:%d,%d:%d] w,h = %d,%d' % (
                x_min,x_max,y_min,y_max,x_max-x_min+1,y_max-y_min+1)
            print ' shift: (%.6f,%.6f) arcsec relative to stamp center' % (
                model.centroid().x,model.centroid().y)

        return cropped_stamp.array[np.newaxis,:,:],cropped_bounds

    @staticmethod
    def add_args(parser):
        """Add command-line arguments for constructing a new :class:`Engine`.

        The added arguments are our constructor parameters with '_' replaced by '-' in the names.
        Note that constructor parameter defaults are specified here rather than in the constructor,
        so that they are included in command-line help.

        Args:
            parser(argparse.ArgumentParser): Arguments will be added to this parser object using its
                add_argument method.
        """
        parser.add_argument('--min-snr', type = float, default = 0.5, metavar = 'SNR',
            help = 'Simulate signals from individual sources down to this S/N threshold.')
        parser.add_argument('--truncate-radius', type = float, default = 30., metavar = 'SIZE',
            help = 'All extended sources are truncated at this radius in arcseconds.')
        parser.add_argument('--no-margin', action = 'store_true',
            help = 'Do not simulate the tails of objects just outside the field.')
        parser.add_argument('--verbose-render', action = 'store_true',
            help = 'Provide verbose output on rendering process.')

    @classmethod
    def from_args(cls,survey,args):
        """Create a new :class:`Engine` object from a set of arguments.

        Args:
            survey(descwl.survey.Survey): Survey that rendered images will simulate.
            args(object): A set of arguments accessed as a :py:class:`dict` using the
                built-in :py:func:`vars` function. Any extra arguments beyond those defined
                in :func:`add_args` will be silently ignored.

        Returns:
            :class:`Engine`: A newly constructed Engine object.
        """
        # Look up the named constructor parameters.
        pnames = (inspect.getargspec(cls.__init__)).args[1:]
        # Get a dictionary of the arguments provided.
        args_dict = vars(args)
        # Filter the dictionary to only include constructor parameters.
        filtered_dict = { key:args_dict[key] for key in (set(pnames) & set(args_dict)) }
        return cls(survey,**filtered_dict)
