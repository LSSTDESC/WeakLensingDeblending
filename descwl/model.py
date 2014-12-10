"""Model astronomical sources.
"""

import inspect

class Galaxy(object):
    """Source model for a galaxy.

    Args:
        disk_flux(float): Total flux in detected electrons of Sersic n=1 component.
        bulge_flux(float): Total flux in detected electrons of Sersic n=4 component.
        agn_flux(float): Total flux in detected electrons of PSF-like component.
    """
    def __init__(self,disk_flux,bulge_flux,agn_flux):
        print 'Galaxy: fluxes =',(disk_flux,bulge_flux,agn_flux)

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
            :class:`Galaxy`: A newly created galaxy source model or None if this object
                should not be simulated.

        Raises:
            RuntimeError: Catalog is missing AB flux value in requested filter band.
        """
        # Calculate the object's total flux in detected electrons.
        try:
            ab_magnitude = entry[filter_band + '_ab']
        except KeyError:
            raise RuntimeError('Catalog entry is missing %s-band AB flux value.')
        total_flux = self.survey.get_flux(ab_magnitude)
        # Calculate the flux of each component in detected electrons.
        total_fluxnorm = entry['fluxnorm_disk'] + entry['fluxnorm_bulge'] + entry['fluxnorm_agn']
        disk_flux = 0. if self.no_disk else entry['fluxnorm_disk']/total_fluxnorm*total_flux
        bulge_flux = 0. if self.no_bulge else entry['fluxnorm_bulge']/total_fluxnorm*total_flux
        agn_flux = 0. if self.no_agn else entry['fluxnorm_agn']/total_fluxnorm*total_flux
        # Is there any flux to simulate?
        if disk_flux + bulge_flux + agn_flux == 0:
            return None
        return Galaxy(disk_flux,bulge_flux,agn_flux)

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
