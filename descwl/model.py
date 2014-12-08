"""Model astronomical sources.
"""

class Galaxy(object):
    """Source model for a galaxy.

    Args:
        total_flux(float): Total source flux 
    """
    def __init__(self,total_flux):
        pass

    @staticmethod
    def from_catalog(entry,args):
        """
        Create a new Galaxy source model from catalog parameters.

        Args:
            entry(ndarray): A numpy array with named columns that specifies the galaxy source
                parameter values to use.
            args(object): An object with named attributes that control how the parameter
                values are used to construct the source model.

        Returns:
            Galaxy: A newly created Galaxy object.
        """
        return Galaxy(total_flux = 1.0)

def add_galaxy_args(parser):
    """Add command-line arguments for building galaxy source models from catalog parameters.

    Args:
        parser(argparse.ArgumentParser): Arguments will be added to this parser object using its
            add_argument method.
    """
    parser.add_argument('--no-disk', action = 'store_true',
        help = 'Ingore any disk (Sersic n=1) component of galaxy sources.')
    parser.add_argument('--no-bulge', action = 'store_true',
        help = 'Ingore any bulge (Sersic n=4) component of galaxy sources.')
