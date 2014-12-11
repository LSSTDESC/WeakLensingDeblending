"""Load source parameters from catalog files.
"""

import math
import inspect

from astropy.io import ascii

class Reader(object):
    """Read source parameters from a catalog to simulate an image.

    The entire catalog is read when the object is constructed, independently of the
    only_id and skip_id parameter values. The recommended method for iterating over
    catalog entries is::

        catalog = descwl.catalog.Reader(...)
        for entry in catalog:
            ...

    Note that constructor parameter defaults are specified in the add_args() function.

    Args:
        catalog_name(str): Name of catalog file, which must exist and be readable.
        ra_center(float): Right ascension of image center in degrees.
        dec_center(float): Declination of image center in degrees.
        only_id(array): Only read ids in this array of integer ids.
        skip_id(array): Skip ids in this array of integer ids.

    Raises:
        RuntimeError: Missing required catalog_name arg.
    """
    def __init__(self,catalog_name,ra_center,dec_center,only_id,skip_id):
        if not catalog_name:
            raise RuntimeError('Missing required catalog_name arg.')
        self.catalog_name = catalog_name
        self.ra_center = ra_center
        self.dec_center = dec_center
        self.only_id = only_id
        self.skip_id = skip_id
        self.table = ascii.read(catalog_name, format='basic')

    def visible_entries(self,survey,render_options):
        """Iterate over visible catalog entries.

        Visibility is determined by the combined survey parameters and rendering options.
        If only_id has any entries, then only the specified ids will be considered. Otherwise,
        all ids are considered. Any ids in skip_id will be silently ignored.

        Args:
            survey(:class:`descwl.survey.Survey`): Survey parameters used to determine which
                entries are visible.
            render_options(:class:`descwl.render.Options`): Rendering options used to determine
                which entries are visible.
        """
        # Calculate the margin size in arcsecs.
        margin_size = 0. if render_options.no_margin else 0.5*render_options.truncate_size
        # Calculate the RA,DEC limits of visible entries in degrees.
        arcsec2deg = 1./3600.
        ra_scale = math.cos(math.radians(self.ra_center))
        ra_size = (0.5*survey.image_width*survey.pixel_scale + margin_size)*arcsec2deg/ra_scale
        ra_min = self.ra_center - 0.5*ra_size
        ra_max = self.ra_center + 0.5*ra_size
        dec_size = (0.5*survey.image_height*survey.pixel_scale + margin_size)*arcsec2deg
        dec_min = self.dec_center - 0.5*dec_size
        dec_max = self.dec_center + 0.5*dec_size
        # Iterate over all catalog entries.
        for row in self.table:
            if self.only_id and row['id'] not in self.only_id:
                continue
            if self.skip_id and row['id'] in self.skip_id:
                continue
            ra,dec = row['ra'],row['dec']
            if ra < ra_min or ra > ra_max or dec < dec_min or dec > dec_max:
                continue
            # If we get this far, the entry is visible.
            yield row

    @staticmethod
    def add_args(parser):
        """Add command-line arguments for constructing a new :class:`Reader`.

        The added arguments are our constructor parameters with '_' replaced by '-' in the names.
        Note that constructor parameter defaults are specified here rather than in the constructor,
        so that they are included in command-line help.

        Args:
            parser(argparse.ArgumentParser): Arguments will be added to this parser object using its
                add_argument method.
        """
        parser.add_argument('--catalog-name', type = str, default = None, metavar = 'NAME',
            help = 'Name of catalog file, which must exist and be readable.')
        parser.add_argument('--ra-center', type = float, default = 0.5, metavar = 'RA',
            help = 'Right ascension of image center in degrees.')
        parser.add_argument('--dec-center', type = float, default = 0.0, metavar = 'DEC',
            help = 'Declination of image center in degrees.')
        parser.add_argument('--only-id', type = int, action = 'append', metavar = 'ID',
            help = 'Use row with ID from the input catalog. May be repeated. All other rows will be ignored.')
        parser.add_argument('--skip-id', type = int, action = 'append', metavar = 'ID',
            help = 'Skip row with ID from the input catalog. May be repeated.')

    @classmethod
    def from_args(cls,args):
        """Create a new :class:`Reader` object from a set of arguments.

        Args:
            args(object): A set of arguments accessed as a :py:class:`dict` using the
                built-in :py:func:`vars` function. Any extra arguments beyond those defined
                in :func:`add_args` will be silently ignored.

        Returns:
            :class:`Reader`: A newly constructed Reader object.
        """
        # Look up the named constructor parameters.
        pnames = (inspect.getargspec(cls.__init__)).args[1:]
        # Get a dictionary of the arguments provided.
        args_dict = vars(args)
        # Filter the dictionary to only include constructor parameters.
        filtered_dict = { key:args_dict[key] for key in (set(pnames) & set(args_dict)) }
        return cls(**filtered_dict)
