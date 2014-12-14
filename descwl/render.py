"""Render source models as simulated survey observations.
"""

import math
import inspect

import numpy as np

import galsim

class Engine(object):
    """Rendering engine to simulate survey observations.

    Args:
        survey(descwl.survey.Survey): Survey that rendered images will simulate.
        min_snr(float): Simulate signals from individual sources down to this S/N threshold,
            where the signal N is calculated for the full exposure time and the noise N is
            set by the expected fluctuations in the sky background during a full exposure.
        truncate_size(float): Truncate sources to a square mask with this full size in arcseconds.
        no_margin(bool): Do not simulate the tails of objects just outside the field.
        verbose_render(bool): Provide verbose output on rendering process.
    """
    def __init__(self,survey,min_snr,truncate_size,no_margin,verbose_render):
        self.survey = survey
        self.min_snr = min_snr
        self.truncate_size = truncate_size
        self.no_margin = no_margin
        self.verbose_render = verbose_render
        # Calculate pixel flux threshold in electrons per pixel that determines how big a
        # bounding box we simulate for each source.
        sky_noise = math.sqrt(survey.get_sky_level())
        self.pixel_cut = self.min_snr*sky_noise
        # Initialize our GalSim parameters.
        self.galsim_params = galsim.GSParams(maximum_fft_size=32768)

    def render_galaxy(self,galaxy):
        """Render a galaxy model for a simulated survey.

        Args:
            galaxy(descwl.model.Galaxy): Model of the galaxy to render.

        Returns:
            :class:`numpy.ndarray`: Array of shape (nstamp,width,height) pixel values that represent
                nstamp postage-stamp images with the same dimensions (width,height) calculated
                based on the rendering options provided. Returns None if this galaxy has no
                pixels above threshold that are visible in the simulated survey image.
        """
        # Calculate the offset of the source center from the bottom-left corner of the
        # simulated image in floating-point pixel units.
        centroid = galaxy.model.centroid()
        x_center_pixels,y_center_pixels = self.survey.get_image_coordinates(centroid.x,centroid.y)

        # Calculate the bounding box extents to use in floating-point pixel units.
        # Use the maximum-sized stamp for now.
        half_width_pixels = 0.5*self.truncate_size/self.survey.pixel_scale
        half_height_pixels = 0.5*self.truncate_size/self.survey.pixel_scale

        # Calculate the bounding box to use for simulating this galaxy.
        x_min = int(math.floor(x_center_pixels - half_width_pixels))
        y_min = int(math.floor(y_center_pixels - half_height_pixels))
        # Use floor instead of ceil here since upper bounds are inclusive.
        x_max = int(math.floor(x_center_pixels + half_width_pixels))
        y_max = int(math.floor(y_center_pixels + half_height_pixels))
        bounds = galsim.BoundsI(x_min,x_max,y_min,y_max)

        # Is there any overlap with the simulated image?
        survey_overlap = bounds & self.survey.image.bounds
        if survey_overlap.area() == 0:
            return None

        # Calculate the offset of the bounding box center from the image center in arcsecs.
        dx_stamp_arcsec = 0.5*(x_min + x_max+1 - self.survey.image_width)*self.survey.pixel_scale
        dy_stamp_arcsec = 0.5*(y_min + y_max+1 - self.survey.image_height)*self.survey.pixel_scale

        # Shift the model to the bounding box center and convolve with the survey PSF.
        model = galsim.Convolve([
            galaxy.model.shift(dx=-dx_stamp_arcsec,dy=-dy_stamp_arcsec),
            self.survey.psf_model
            ],gsparams=self.galsim_params)

        # Render the model in its own postage stamp.
        stamp = galsim.Image(bounds = bounds,scale = self.survey.pixel_scale, dtype = np.float64)
        model.drawImage(image = stamp, use_true_center = True)        
        self.survey.image[survey_overlap] += stamp[survey_overlap]

        # Draw directly into the survey image, with no truncation.
        #model.drawImage(image = self.survey.image,add_to_image = True,use_true_center = True)

        if self.verbose_render:
            print 'Rendering galaxy model for id = %d with z = %.3f' % (
                galaxy.identifier,galaxy.redshift)
            print 'bounds: [%d:%d,%d:%d] w,h = %d,%d' % (
                x_min,x_max,y_min,y_max,x_max-x_min+1,y_max-y_min+1)
            print ' shift: (%.6f,%.6f) arcsec relative to stamp center' % (
                model.centroid().x,model.centroid().y)

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
        parser.add_argument('--truncate-size', type = float, default = 30., metavar = 'SIZE',
            help = 'Truncate sources to a square mask with this full size in arcseconds.')
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
