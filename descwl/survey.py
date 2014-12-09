"""Manage the parameters that define a simulated survey's camera design and observing conditions.
"""

class Survey(object):
    """Survey camera and observing parameters.

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
    """
    def __init__(self,image_width,image_height,pixel_scale,exposure_time,zero_point,instrumental_psf_fwhm,
        zenith_psf_fwhm,atmospheric_psf_beta,sky_brightness,airmass,extinction):
        self.image_width = image_width
        self.image_height = image_height
        self.pixel_scale = pixel_scale
        self.exposure_time = exposure_time
        self.zero_point = zero_point
        self.instrumental_psf_fwhm = instrumental_psf_fwhm
        self.zenith_psf_fwhm = zenith_psf_fwhm
        self.atmospheric_psf_beta = atmospheric_psf_beta
        self.sky_brightness = sky_brightness
        self.airmass = airmass
        self.extinction = extinction

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
    def add_args(parser):
        """Add command-line arguments for constructing a new Survey.

        The added arguments are our constructor parameters with '_' replaced by '-' in the names.
        Note that constructor parameter defaults are specified here rather than in the constructor,
        so that they are included in command-line help.

        Args:
            parser(argparse.ArgumentParser): Arguments will be added to this parser object using its
                add_argument method.
        """
        parser.add_argument('--survey', choices = ['LSST','DES','CFHT'], default = 'LSST',
            help = 'Use default camera and observing parameters appropriate for this survey.')
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

    @classmethod
    def from_args(cls,args):
        """Create a new Survey object from a set of arguments.

        Args:
            args(args): A set of arguments convertible to a dictionary via the built-in vars() method.
                The argument set should include those defined in add_args(). A filter_band argument
                specifying the LSST band to use for setting defaults is also required, taking one of
                the values 'u','g','r','i','z','y'. Not all (survey,filter) combinations are supported.

        Returns:
            Survey: A newly constructed Survey object.

        Raises:
            RuntimeError: defaults not available for the requested survey and filter band.
        """
        args_dict = vars(args)
        survey = args_dict['survey']
        filter_band = args_dict['filter_band']
        # Do we have defaults for the requested survey and filter band?
        if survey not in Survey._defaults:
            raise RuntimeError('Defaults not defined yet for %s.' % survey)
        survey_defaults = Survey._defaults[survey]
        if filter_band not in survey_defaults:
            raise RuntimeError('Defaults not defined yet for %s %s-band.' % (survey,filter_band))
        # Set defaults.
        ctor_params = Survey._defaults['*']
        if '*' in survey_defaults:
            ctor_params.update(survey_defaults['*'])
        ctor_params.update(survey_defaults[filter_band])
        # Override defaults using the args provided.
        for name in ctor_params:
            if name in args_dict and args_dict[name] is not None:
                ctor_params[name] = args_dict[name]
        return Survey(**ctor_params)
