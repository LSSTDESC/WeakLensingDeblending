"""Render source models as simulated survey observations.
"""

import inspect

class Options(object):
    """Options to control rendering of sources as simulated survey observations.

    Args:
        min_snr(float): Simulate signals from individual sources down to this S/N threshold,
            where the signal N is calculated for the full exposure time and the noise N is
            set by the expected fluctuations in the sky background during a full exposure.
    """
    def __init__(self,min_snr):
        pass

    @staticmethod
    def add_args(parser):
        """Add command-line arguments for constructing a new :class:`Options`.

        The added arguments are our constructor parameters with '_' replaced by '-' in the names.
        Note that constructor parameter defaults are specified here rather than in the constructor,
        so that they are included in command-line help.

        Args:
            parser(argparse.ArgumentParser): Arguments will be added to this parser object using its
                add_argument method.
        """
        parser.add_argument('--min-snr', type = float, default = 0.5, metavar = 'SNR',
            help = 'Simulate signals from individual sources down to this S/N threshold.')

    @classmethod
    def from_args(cls,args):
        """Create a new :class:`Options` object from a set of arguments.

        Args:
            args(object): A set of arguments accessed as a :py:class:`dict` using the
                built-in :py:func:`vars` function. Any extra arguments beyond those defined
                in :func:`add_args` will be silently ignored.

        Returns:
            :class:`Options`: A newly constructed Options object.
        """
        # Look up the named constructor parameters.
        pnames = (inspect.getargspec(cls.__init__)).args[1:]
        # Get a dictionary of the arguments provided.
        args_dict = vars(args)
        # Filter the dictionary to only include constructor parameters.
        filtered_dict = { key:args_dict[key] for key in (set(pnames) & set(args_dict)) }
        return cls(**filtered_dict)
