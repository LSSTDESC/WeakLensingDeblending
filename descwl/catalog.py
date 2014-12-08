"""Load source parameters from catalog files.
"""

import inspect

from astropy.io import ascii

class Reader(object):
    """Read source parameters from a catalog to simulate an image.

    The entire catalog is read when the object is constructed, independently of the
    only_row and skip_row parameter values. The recommended method for iterating over
    catalog entries is::

        catalog = descwl.catalog.Reader(...)
        for entry in catalog:
            ...

    Note that constructor parameter defaults are specified in the add_args() function.

    Args:
        catalog_name(str): Name of catalog file, which must exist and be readable.
        ra_center(float): Right ascension of image center in degrees.
        dec_center(float): Declination of image center in degrees.
        only_row(array): Only read rows in this array of integer row indices.
        skip_row(array): Skip rows in this array of integer row indices.
    """
    def __init__(self,catalog_name,ra_center,dec_center,only_row,skip_row):
        self.only_row = only_row
        self.skip_row = skip_row
        self.table = ascii.read(catalog_name, format='basic')

    def __iter__(self):
        """Iterate over catalog entries.

        If only_row has any entries, then only the specified rows will be considered. Otherwise,
        all rows are considered. Any rows in skip_row will be silently ignored.
        """
        for index,row in enumerate(self.table):
            if self.only_row and index not in self.only_row:
                continue
            if self.skip_row and index in self.skip_row:
                continue
            yield row

    @staticmethod
    def add_args(parser):
        """Add command-line arguments for constructing a new CatalogReader.

        The added arguments are our constructor parameters. Note that constructor parameter defaults
        are specified here rather than in the constructor, so that they are included in command-line
        help.

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
        parser.add_argument('--only-row', type = int, action = 'append', metavar = 'N',
            help = 'Use row N from the input catalog. May be repeated. All other lines will be ignored.')
        parser.add_argument('--skip-row', type = int, action = 'append', metavar = 'N',
            help = 'Skip row N from the input catalog. May be repeated.')

    @staticmethod
    def from_args(args):
        """Create a new Reader object from a set of arguments.

        Args:
            args(args): A set of arguments convertible to a dictionary via the built-in vars() method.

        Returns:
            Reader: A newly constructed Reader object.
        """
        # Look up the named constructor parameters.
        pnames = (inspect.getargspec(Reader.__init__)).args[1:]
        # Get a dictionary of the arguments provided.
        args_dict = vars(args)
        # Filter the dictionary to only include constructor parameters.
        filtered_dict = { key:args_dict[key] for key in (set(pnames) & set(args_dict)) }
        print filtered_dict
        return Reader(**filtered_dict)
