"""Manage the parameters that define a simulated survey's camera design and observing conditions.
"""

import math

import galsim

class Survey(object):
    """Survey camera and observing parameters.

        No default values are assigned to constructor args since these are handled at run time
        based on a requested (survey_name,filter_band) combination using :func:`get_defaults`.

    Args:
        survey_name(str): Use default camera and observing parameters appropriate for this survey.
        filter_band(str): LSST imaging band to simulate. Must be one of 'u','g','r','i','z','y'.
        image_width(int): Simulated camera image width in pixels.
        image_height(int): Simulated camera image height in pixels.
        pixel_scale(float): Simulated camera pixel scale in arcseconds per pixel.
        exposure_time(float): Simulated camera total exposure time seconds.
        zero_point(float): Simulated camera zero point in electrons per second at 24th magnitude.
        mirror_diameter(float): Size of the primary mirror in meters to use for the optical PSF,
            or zero if no optical PSF should be simulated.
        obscuration_fraction(float): Linear obscuration fraction of the primary mirror to
            use for the optical PSF. Ignored if mirror_diameter is zero.
        zenith_psf_fwhm(float): FWHM of the atmospheric PSF at zenith in arcseconds.
        atmospheric_psf_beta(float): Moffat beta parameter of the atmospheric PSF, or use a Kolmogorov
            PSF if beta <= 0.
        sky_brightness(float): Sky brightness in mags/sq.arcsec during the observation.
        airmass(float): Optical path length through the atmosphere relative to the zenith path length.
        extinction(float): Exponential exctinction coefficient for atmospheric absorption.

    Raises:
        RuntimeError: Missing or extra arguments provided.
    """
    def __init__(self,**args):
        if set(args.keys()) != set(Survey._parameter_names):
            raise RuntimeError('Missing or extra arguments provided to Survey constructor.')
        self.args = args
        self.__dict__.update(args)
        # Build our atmospheric PSF model.
        atmospheric_psf_fwhm = self.zenith_psf_fwhm*self.airmass**0.6
        if self.atmospheric_psf_beta > 0:
            atmospheric_psf_model = galsim.Moffat(
                beta = self.atmospheric_psf_beta, fwhm = atmospheric_psf_fwhm)
        else:
            atmospheric_psf_model = galsim.Kolmogorov(fwhm = atmospheric_psf_fwhm)
        # Combine with our optical PSF model, if any.
        if self.mirror_diameter > 0:
            lambda_over_diameter = 3600*math.degrees(
                1e-10*Survey._central_wavelength[self.filter_band]/self.mirror_diameter)
            assert obscuration >= 0. and obscuration <= 1., 'Invalid obscuration_fraction'
            optical_psf_model = galsim.Airy(lam_over_diam = lambda_over_diameter,
                obscuration = self.obscuration_fraction)
            self.psf_model = galsim.Add(atmospheric_psf_model,optical_psf_model)
        else:
            self.psf_model = atmospheric_psf_model
        # Create an empty image using (0,0) to index the lower-left corner pixel.
        self.image_bounds = galsim.BoundsI(0,self.image_width-1,0,self.image_height-1)
        self.image = galsim.Image(bounds = self.image_bounds,scale=self.pixel_scale,
            dtype = np.float64)

    def description(self):
        """Describe the survey we simulate.

        Returns:
            str: Description of the camera design and observing conditions we simulate.
        """
        return 'Simulating %s %s-band survey with %r' % (
            self.survey_name,self.filter_band,self.args)

    def get_flux(self,ab_magnitude):
        """Convert source magnitude to flux.

        The calculation includes the effects of atmospheric extinction.

        Args:
            ab_magnitude(float): AB magnitude of source.

        Returns:
            float: Flux in detected electrons.
        """
        ab_magnitude += self.extinction*(self.airmass - 1.)
        return self.exposure_time*self.zero_point*10**(-0.4*(ab_magnitude-24))

    def get_sky_level(self):
        """Calculate the mean sky background level.

        Returns:
            float: Mean sky background flux in detected electrons per pixel.
        """
        return self.get_flux(self.sky_brightness)*self.pixel_scale**2

    def get_image_coordinates(self,dx_arcsecs,dy_arcsecs):
        """Convert a physical offset from the image center into image coordinates.

        Args:
            dx_arcsecs(float): Offset from the image center in arcseconds.
            dy_arcsecs(float): Offset from the image center in arcseconds.

        Returns:
            tuple: Corresponding floating-point image coordinates (x_pixels,y_pixels)
                whose :func:`math.floor` value gives pixel indices and whose...
        """
        x_pixels = 0.5*self.image_width + dx_arcsecs/self.pixel_scale
        y_pixels = 0.5*self.image_height + dy_arcsecs/self.pixel_scale
        return x_pixels,y_pixels

    # Survey constructor parameter names. The order established here is used by print_defaults().
    _parameter_names = (
        'survey_name','filter_band',
        'image_width','image_height','pixel_scale','exposure_time','zero_point',
        'mirror_diameter','obscuration_fraction',
        'zenith_psf_fwhm','atmospheric_psf_beta','sky_brightness',
        'airmass','extinction'
        )

    # Central wavelengths in Angstroms for each LSST filter band, calculated from the
    # baseline total filter throughputs tabulated at
    # http://dev.lsstcorp.org/cgit/LSST/sims/throughputs.git/snapshot/throughputs-1.2.tar.gz
    _central_wavelength = {
        'u':3592.13, 'g':4789.98, 'r':6199.52, 'i':7528.51, 'z':8689.83, 'y':9674.05
        }

    # Default constructor arg values for different (survey,filter_band) combinations.
    _defaults = {
        '*': {
            'mirror_diameter': 0.0,
            'obscuration_fraction': 0.0,
            'atmospheric_psf_beta': 0.0,
            'airmass': 1.2,
        },
        'LSST': {
            '*': {
                'image_width': 4096,
                'image_height': 4096,
                'pixel_scale': 0.2,
            },
            'i': {
                'exposure_time': 6900.,
                'sky_brightness': 20.0,
                'zenith_psf_fwhm': 0.67,
                'zero_point': 41.5,
                'extinction': 0.07,
            },
            'r': {
                'exposure_time': 6900.,
                'sky_brightness': 21.3,
                'zenith_psf_fwhm': 0.70,
                'zero_point': 55.8,
                'extinction': 0.10,
            },
        },
        'DES': {
            '*': {
                'image_width': 3115,
                'image_height': 3115,
                'pixel_scale': 0.263,
            },
            'i': {
                'exposure_time': 1000.,
                'sky_brightness': 20.1,
                'zenith_psf_fwhm': 0.79,
                'zero_point': 12.5,
                'extinction': 0.07,
            },
            'r' : {
                'exposure_time': 800.,
                'sky_brightness': 21.1,
                'zenith_psf_fwhm': 0.79,
                'zero_point': 16.8,
                'extinction': 0.10,
            },
        },
        'CFHT': {
            '*': {
                'image_width': 4428,
                'image_height': 4428,
                'pixel_scale': 0.185,
            },
            'i': {
                'exposure_time': 4300.,
                'sky_brightness': 20.3,
                'zenith_psf_fwhm': 0.64,
                'zero_point': 10.0,
                'extinction': 0.07,
            },
            'r' : {
                'exposure_time': 2000.,
                'sky_brightness': 20.8,
                'zenith_psf_fwhm': 0.71,
                'zero_point': 13.5,
                'extinction': 0.10,
            },
        }
    }

    @staticmethod
    def print_defaults():
        """Print parameters for all available (survey,filter band) combinations.
        """
        for survey_name,survey_defaults in Survey._defaults.iteritems():
            if survey_name == '*':
                continue
            for filter_band,combo_defaults in survey_defaults.iteritems():
                if filter_band == '*':
                    continue
                defaults = Survey.get_defaults(survey_name,filter_band)
                print '%s %s-band: %r' % (survey_name,filter_band,defaults)

    @staticmethod
    def add_args(parser):
        """Add command-line arguments for constructing a new :class:`Survey`.

        The added arguments are our constructor parameters with '_' replaced by '-' in the names.
        No default values are assigned to the added args since these are handled at run time
        based on a requested (survey_name,filter_band) combination using :func:`get_defaults`.

        Args:
            parser(argparse.ArgumentParser): Arguments will be added to this parser object using its
                add_argument method.
        """
        parser.add_argument('--survey-name', choices = ['LSST','DES','CFHT'], default = 'LSST',
            help = 'Use default camera and observing parameters appropriate for this survey.')
        parser.add_argument('--filter-band', choices = ['u','g','r','i','z','y'], default = 'i',
            help = 'LSST imaging band to simulate')
        parser.add_argument('--image-width', type = int, metavar = 'W',
            help = 'Simulated camera image width in pixels.')
        parser.add_argument('--image-height', type = int, metavar = 'H',
            help = 'Simulated camera image height in pixels.')
        parser.add_argument('--pixel-scale', type = float, metavar = 'S',
            help = 'Simulated camera pixel scale in arcseconds per pixel.')
        parser.add_argument('--exposure-time', type = float, metavar = 'T',
            help = 'Simulated camera total exposure time seconds.')
        parser.add_argument('--zero-point', type = float, metavar = 's0',
            help = 'Simulated camera zero point in electrons per second at 24th magnitude.')
        parser.add_argument('--mirror-diameter', type = float, metavar = 'D',
            help = 'Size of the primary mirror in meters for the optical PSF (or zero for no PSF).')
        parser.add_argument('--obscuration-fraction', type = float, metavar = 'F',
            help = 'Linear obscuration fraction of the primary mirror for the optical PSF.')
        parser.add_argument('--zenith-psf-fwhm', type = float, metavar = 'FWHM',
            help = 'FWHM of the atmospheric PSF at zenith in arcseconds.')
        parser.add_argument('--atmospheric-psf-beta', type = float, metavar = 'BETA',
            help = 'Moffat beta parameter of the atmospheric PSF, or use a Kolmogorov PSF if beta <= 0.')
        parser.add_argument('--sky-brightness', type = float, metavar = 'B',
            help = 'Sky brightness in mags/sq.arcsec during the observation.')
        parser.add_argument('--airmass', type = float, metavar = 'X',
            help = 'Optical path length through the atmosphere relative to the zenith path length.')
        parser.add_argument('--extinction', type = float, metavar = 'k',
            help = 'Exponential exctinction coefficient for atmospheric absorption.')

    @staticmethod
    def get_defaults(survey_name,filter_band):
        """Get survey camera and observing parameter default values.

        Args:
            survey_name(str): Use defaults for this survey. Case is significant.
            filter_band(str): Use defaults for this filter band. Case is significant.

        Returns:
            dict: Dictionary of parameter (name,value) pairs.

        Raises:
            RuntimeError: Defaults not yet defined for requested combination.
        """
        # Do we have defaults for the requested survey name and filter band?
        if survey_name not in Survey._defaults:
            raise RuntimeError('Defaults not defined yet for %s.' % survey_name)
        survey_defaults = Survey._defaults[survey_name]
        if filter_band not in survey_defaults:
            raise RuntimeError('Defaults not defined yet for %s %s-band.' % (survey_name,filter_band))
        # Set defaults.
        defaults = Survey._defaults['*']
        if '*' in survey_defaults:
            defaults.update(survey_defaults['*'])
        defaults.update(survey_defaults[filter_band])
        return defaults

    @classmethod
    def from_args(cls,args):
        """Create a new :class:`Survey` object from a set of arguments.

        Args:
            args(object): A set of arguments accessed as a :py:class:`dict` using the
                built-in :py:func:`vars` function. The argument set must include those defined in
                :func:`add_args`. Two additional arguments are also required: survey_name and
                filter_band. These establish the constructor parameter defaults via :func:get_defaults.
                Any extra arguments beyond those defined in :func:`add_args`, survey_name, and
                filter_band will be silently ignored.

        Returns:
            :class:`Survey`: A newly constructed Survey object.

        Raises:
            RuntimeError: Defaults not yet defined for requested combination.
        """
        args_dict = vars(args)
        survey_name = args_dict['survey_name']
        filter_band = args_dict['filter_band']
        # Set defaults for this (survey,band) combination.
        ctor_params = Survey.get_defaults(survey_name,filter_band)
        # Override defaults using the args provided.
        for parameter_name in ctor_params:
            if parameter_name in args_dict and args_dict[parameter_name] is not None:
                ctor_params[parameter_name] = args_dict[parameter_name]
        return Survey(survey_name=survey_name,filter_band=filter_band,**ctor_params)
