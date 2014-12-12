"""Model astronomical sources.
"""

import math
import inspect

class Galaxy(object):
    """Source model for a galaxy.

    Galaxies are modeled using up to three components: a disk (Sersic n=1), a bulge
    (Sersic n=4), and an AGN (PSF-like). Not all components are required.  All components
    are assumed to have the same centroid and the extended (disk,bulge) components are
    assumed to have the same position angle.

    Args:
        dx_arcsecs(float): Horizontal offset of catalog entry's centroid from image center
            in arcseconds.
        dy_arcsecs(float): Vertical offset of catalog entry's centroid from image center
            in arcseconds.
        beta_radians(float): Position angle beta of Sersic components in radians, measured
            anti-clockwise from the positive x-axis. Ignored if disk_flux and bulge_flux are
            both zero.
        disk_flux(float): Total flux in detected electrons of Sersic n=1 component.
        disk_hlr_arcsecs(float): Half-light radius sqrt(a*b) of circularized 50% isophote
            for Sersic n=1 component, in arcseconds. Ignored if disk_flux is zero.
        disk_q(float): Ratio b/a of 50% isophote semi-minor (b) to semi-major (a) axis
            lengths for Sersic n=1 component. Ignored if disk_flux is zero.
        bulge_flux(float): Total flux in detected electrons of Sersic n=4 component.
        bulge_hlr_arcsecs(float): Half-light radius sqrt(a*b) of circularized 50% isophote
            for Sersic n=4 component, in arcseconds. Ignored if bulge_flux is zero.
        bulge_q(float): Ratio b/a of 50% isophote semi-minor (b) to semi-major (a) axis
            lengths for Sersic n=4 component. Ignored if bulge_flux is zero.
        agn_flux(float): Total flux in detected electrons of PSF-like component.
    """
    def __init__(self,dx_arcsecs,dy_arcsecs,beta_radians,disk_flux,disk_hlr_arcsecs,disk_q,
            bulge_flux,bulge_hlr_arcsecs,bulge_q,agn_flux):
        print 'Galaxy:',(disk_flux,bulge_flux,agn_flux)

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

    def from_catalog(self,entry,dx_arcsecs,dy_arcsecs,filter_band):
        """Build a :class:Galaxy object from a catalog entry.

        Fluxes are distributed between the three possible components (disk,bulge,AGN) assuming
        that each component has the same spectral energy distribution, so that the resulting
        proportions are independent of the filter band.

        Args:
            entry(astropy.table.Row): A single row from a galaxy :mod:`descwl.catalog`.
            dx_arcsecs(float): Horizontal offset of catalog entry's centroid from image center
                in arcseconds.
            dy_arcsecs(float): Vertical offset of catalog entry's centroid from image center
                in arcseconds.
            filter_band(str): The LSST filter band to use for calculating flux, which must
                be one of 'u','g','r','i','z','y'.

        Returns:
            :class:`Galaxy`: A newly created galaxy source model or None if this object
                should not be simulated.

        Raises:
            RuntimeError: Catalog entry is missing AB flux value in requested filter band.
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
        # Calculate the position of angle of the Sersic components, which are assumed to be the same.
        if disk_flux > 0:
            beta_radians = math.radians(entry['pa_disk'])
            if bulge_flux > 0:
                assert (entry['pa_disk'] == entry['pa_bulge'],
                    'Disk, bulge components have different beta.')
        elif bulge_flux > 0:
            beta_radians = math.radians(entry['pa_bulge'])
        else:
            # This might happen if we only have an AGN component.
            beta_radians = None
        # Calculate shapes hlr = sqrt(a*b) and q = b/a of Sersic components.
        if disk_flux > 0:
            a_d,b_d = entry['a_d'],entry['b_d']
            disk_hlr_arcsecs = math.sqrt(a_d*b_d)
            disk_q = b_d/a_d
        else:
            disk_hlr_arcsecs,disk_q = None,None
        if bulge_flux > 0:
            a_b,b_b = entry['a_b'],entry['b_b']
            bulge_hlr_arcsecs = math.sqrt(a_b*b_b)
            bulge_q = b_b/a_b
            bulge_beta = math.radians(entry['pa_bulge'])
        else:
            bulge_hlr_arcsecs,bulge_q = None,None
        return Galaxy(dx_arcsecs,dy_arcsecs,beta_radians,disk_flux,disk_hlr_arcsecs,disk_q,
            bulge_flux,bulge_hlr_arcsecs,bulge_q,agn_flux)

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