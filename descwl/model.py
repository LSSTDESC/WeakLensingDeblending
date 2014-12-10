"""Model astronomical sources.
"""

import inspect

class Galaxy(object):
    """Source model for a galaxy.

    Args:
        total_flux(float): Total source flux in detected electrons.
    """
    def __init__(self,total_flux):
        print 'Galaxy: total_flux =',total_flux,'elec'

class GalaxyBuilder(object):
    """Build galaxy source models.

    Args:
        survey(descwl.survey.Survey): Survey parameters to use for flux normalization.
        no_disk(bool): Ignore any Sersic n=1 component in the model if it is present in the catalog.
        no_bulge(bool): Ignore any Sersic n=4 component in the model if it is present in the catalog.
        no_agn(bool): Ignore any PSF-like component in the model if it is present in the catalog.
    """
    def __init__(self,survey,no_disk,no_bulge,no_agn):
        if no_disk and no_bulge and no_agn:
            raise RuntimeError('Must build at least one galaxy component.')
        self.survey = survey
        self.no_disk = no_disk
        self.no_bulge = no_bulge
        self.no_agn = no_agn

    def from_catalog(self,entry,filter_band):
        """Build a :class:Galaxy object from a catalog entry.

        Args:
            entry(astropy.table.Row): A single row from a galaxy :mod:`descwl.catalog`.

        Returns:
            :class:`Galaxy`: A newly created galaxy source model.

        Raises:
            RuntimeError: Catalog is missing AB flux value in requested filter band.
        """
        try:
            ab_magnitude = entry[filter_band + '_ab']
        except KeyError:
            raise RuntimeError('Catalog entry is missing %s-band AB flux value.')
        total_flux = self.survey.get_flux(ab_magnitude)
        return Galaxy(total_flux)

    @staticmethod
    def add_args(parser):
        """Add command-line arguments for constructing a new :class:`GalaxyBuilder`.

        The added arguments are our constructor parameters with '_' replaced by '-' in the names.

        Args:
            parser(argparse.ArgumentParser): Arguments will be added to this parser object using its
                add_argument method.
        """
        parser.add_argument('--no-disk', action = 'store_true',
            help = 'Ignore any Sersic n=1 component in the model if it is present in the catalog.')
        parser.add_argument('--no-bulge', action = 'store_true',
            help = 'Ignore any Sersic n=4 component in the model if it is present in the catalog.')
        parser.add_argument('--no-agn', action = 'store_true',
            help = 'Ignore any PSF-like component in the model if it is present in the catalog.')        

    @classmethod
    def from_args(cls,survey,args):
        """Create a new :class:`GalaxyBuilder` object from a set of arguments.

        Args:
            args(object): A set of arguments accessed as a :py:class:`dict` using the
                built-in :py:func:`vars` function. Any extra arguments beyond those defined
                in :func:`add_args` will be silently ignored.

        Returns:
            :class:`GalaxyBuilder`: A newly constructed Reader object.
        """
        # Look up the named constructor parameters.
        pnames = (inspect.getargspec(cls.__init__)).args[1:]
        # Get a dictionary of the arguments provided.
        args_dict = vars(args)
        # Filter the dictionary to only include constructor parameters.
        filtered_dict = { key:args_dict[key] for key in (set(pnames) & set(args_dict)) }
        return cls(survey,**filtered_dict)
