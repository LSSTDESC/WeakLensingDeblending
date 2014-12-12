"""Render source models as simulated survey observations.
"""

import math
import inspect

class Engine(object):
    """Rendering engine to simulate survey observations.

    Args:
        min_snr(float): Simulate signals from individual sources down to this S/N threshold,
            where the signal N is calculated for the full exposure time and the noise N is
            set by the expected fluctuations in the sky background during a full exposure.
        truncate_size(float): Truncate sources to a square mask with this full size in arcseconds.
        no_margin(bool): Do not simulate the tails of objects just outside the field.
    """
    def __init__(self,min_snr,truncate_size,no_margin):
        self.min_snr = min_snr
        self.truncate_size = truncate_size
        self.no_margin = no_margin

    def render_galaxy(self,model,survey):
        """Render a galaxy model for a simulated survey.

        Args:
            model(descwl.model.Galaxy): Model of the galaxy to render.
            survey(deswl.survey.Survey): Survey camera and observing parameters to simulate.

        Returns:
            :class:`numpy.ndarray`: Array of shape (nstamp,width,height) pixel values that represent
                nstamp postage-stamp images with the same dimensions (width,height) calculated
                based on the rendering options provided.
        """
        # Calculate pixel flux threshold in electrons per pixel that determines our bounding box.
        sky_noise = math.sqrt(survey.get_sky_level())
        pixel_cut = self.min_snr*sky_noise
        return None

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
        parser.add_argument('--truncate-size', type = float, default = 20., metavar = 'SIZE',
            help = 'Truncate sources to a square mask with this full size in arcseconds.')
        parser.add_argument('--no-margin', action = 'store_true',
            help = 'Do not simulate the tails of objects just outside the field.')

    @classmethod
    def from_args(cls,args):
        """Create a new :class:`Engine` object from a set of arguments.

        Args:
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
        return cls(**filtered_dict)
