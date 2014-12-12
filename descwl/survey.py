"""Manage the parameters that define a simulated survey's camera design and observing conditions.
"""

class Survey(object):
    """Survey camera and observing parameters.

        No default values are assigned to constructor args since these are handled at run time
        based on a requested (survey_name,filter_band) combination using :func:`get_defaults`.

    Args:
        image_width(int): Simulated camera image width in pixels.
        image_height(int): Simulated camera image height in pixels.
        pixel_scale(float): Simulated camera pixel scale in arcseconds per pixel.
        exposure_time(float): Simulated camera total exposure time seconds.
        zero_point(float): Simulated camera zero point in electrons per second at 24th magnitude.
        instrumental_psf_fwhm(float): FWHM of the simulated camera PSF in arcseconds.
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

    def get_flux(self,ab_magnitude):
        """Convert source magnitude to flux.

        Args:
            ab_magnitude(float): AB magnitude of source.

        Returns:
            float: Flux in detected electrons.
        """
        return self.exposure_time*self.zero_point*10**(-0.4*(ab_magnitude-24))

    def get_sky_level(self):
        """Calculate the mean sky background level.

        Returns:
            float: Mean sky background flux in detected electrons per pixel.
        """
        return self.get_flux(self.sky_brightness)*self.pixel_scale**2

    # Survey constructor parameter names. The order established here is used by print_defaults().
    _parameter_names = (
        'image_width','image_height','pixel_scale','exposure_time','zero_point',
        'instrumental_psf_fwhm','zenith_psf_fwhm','atmospheric_psf_beta','sky_brightness',
        'airmass','extinction'
        )

    # Default constructor arg values for different (survey,filter_band) combinations.
    _defaults = {
        '*': {
            'instrumental_psf_fwhm': 0.4,
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
        parser.add_argument('--instrumental-psf-fwhm', type = float, metavar = 'FWHM',
            help = 'FWHM of the simulated camera PSF in arcseconds.')
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
        return Survey(**ctor_params)
