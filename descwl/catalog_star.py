"""Load source parameters from catalog files.

There is a separate :doc:`catalog page </catalog>` with details on the expected catalog
contents and formatting.
"""

import math
import inspect
import os.path

import astropy.table


class Reader(object):
    """Read source parameters from a catalog to simulate an image.

    The entire catalog is read when the object is constructed, independently of the
    only_id and skip_id parameter values. Note that constructor parameter defaults
    are specified in the add_args() function.

    Details on the catalog contents and format are documented on the
    :doc:`catalog page </catalog>`.

    Args:
        catalog_name(str): Name of catalog file, which must exist and be readable. The
            catalog will be read as a FITS table if catalog_name ends with '.fits'.
            Otherwise, it will be read as an ASCII file.
        ra_center(float): Right ascension of image center in degrees.
        dec_center(float): Declination of image center in degrees.
        only_id(array): Only read ids in this array of integer ids.
        skip_id(array): Skip ids in this array of integer ids.

    Raises:
        RuntimeError: Missing required catalog_name arg.
    """
    def __init__(self,catalog_name,ra_center = 0.0,dec_center = 0.0,only_id = [],skip_id = []):
        if not catalog_name:
            raise RuntimeError('Missing required catalog_name arg.')
        
        self.catalog_name = catalog_name
        self.ra_center = ra_center
        self.dec_center = dec_center
        self.only_id = only_id
        self.skip_id = skip_id
        name,ext = os.path.splitext(catalog_name)
        
        if ext == '.fits':
            self.table = astropy.table.Table.read(catalog_name, format = 'fits')
        else:
            self.table = astropy.table.Table.read(catalog_name, format='ascii.basic')
        
    def potentially_visible_entries(self,survey,render_options):
        """Iterate over potentially visible catalog entries.

        Potential visibility is determined by the combined survey parameters and rendering
        options and implemented so that all actually visible entries are included. Returned
        entries might not be actually visible, for example, if they are too faint.
        If only_id has any entries, then only the specified ids will be considered. Otherwise,
        all ids are considered. Any ids in skip_id will be silently ignored.

        Filtering on (ra,dec) to determine visibility is currently not very robust and assumes
        that ra_center is close to zero and that subtracting 360 from any ra > 180 is sufficient
        for interval tests.

        Args:
            survey(:class:`descwl.survey.Survey`): Survey parameters used to determine which
                entries are visible.
            render_options(:class:`descwl.render.Options`): Rendering options used to determine
                which entries are visible.

        Yields:
            tuple: Tuple (entry,dx,dy) of the current catalog entry (:class:`astropy.table.Row`)
                and the x,y offsets (:class:`float`) of this entry's centroid from the image
                center in pixels.
        """
        # Calculate the margin size in arcsecs.
        margin_size = 0. if render_options.no_margin else render_options.truncate_radius
        # Calculate the RA,DEC limits of visible entries in degrees.
        arcsec2deg = 1./3600.
        # Calculate scaling between RA and vertical angular separations in the image.
        ra_scale = math.cos(math.radians(self.dec_center))
        ra_size = (0.5*survey.image_width*survey.pixel_scale + margin_size)*arcsec2deg/ra_scale
        ra_min = self.ra_center - ra_size
        ra_max = self.ra_center + ra_size
        dec_size = (0.5*survey.image_height*survey.pixel_scale + margin_size)*arcsec2deg
        dec_min = self.dec_center - dec_size
        dec_max = self.dec_center + dec_size
        # Iterate over all catalog entries.
        for entry in self.table:
            if self.only_id and entry['id'] not in self.only_id:
                continue
            if self.skip_id and entry['id'] in self.skip_id:
                continue
            ra,dec = entry['ra'],entry['dec']
            if ra > 180:
                ra -= 360.
            if ra < ra_min or ra > ra_max or dec < dec_min or dec > dec_max:
                continue
            # If we get this far, the entry is visible.
            dx_arcsecs = (ra - self.ra_center)/arcsec2deg*ra_scale
            dy_arcsecs = (dec - self.dec_center)/arcsec2deg
            yield entry,dx_arcsecs,dy_arcsecs
        
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
        parser.add_argument('--ra-center', type = float, default = 0.0, metavar = 'RA',
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


class ReaderStar(object):
    """Read source parameters from a catalog to simulate an image.

    The entire catalog is read when the object is constructed, independently of the
    only_id and skip_id parameter values. Note that constructor parameter defaults
    are specified in the add_args() function.

    Details on the catalog contents and format are documented on the
    :doc:`catalog page </catalog>`.

    Args:
        star_catalog_name(str): Name of catalog file, which must exist and be readable. The
            star catalog will be read as a FITS table if catalog_name ends with '.fits'.
            Otherwise, it will be read as an ASCII file.
        ra_center(float): Right ascension of image center in degrees.
        dec_center(float): Declination of image center in degrees.
        only_id(array): Only read ids in this array of integer ids.
        skip_id(array): Skip ids in this array of integer ids.

    Raises:
        RuntimeError: Missing required catalog_name arg.
    """
    def __init__(self,star_catalog_name,ra_center = 0.0,dec_center = 0.0,only_id = [],skip_id = []):
        if not star_catalog_name:
            raise RuntimeError('Missing required catalog_name arg.')
        
        self.star_catalog_name = star_catalog_name
        self.ra_center = ra_center
        self.dec_center = dec_center
        self.only_id = only_id
        self.skip_id = skip_id
        name,ext = os.path.splitext(star_catalog_name)
        
        if ext == '.fits':
            self.table = astropy.table.Table.read(star_catalog_name, format = 'fits')
        else:
            self.table = astropy.table.Table.read(star_catalog_name, format='ascii.basic')
        
    def potentially_visible_entries(self,survey,render_options):
        """Iterate over potentially visible catalog entries.

        Potential visibility is determined by the combined survey parameters and rendering
        options and implemented so that all actually visible entries are included. Returned
        entries might not be actually visible, for example, if they are too faint.
        If only_id has any entries, then only the specified ids will be considered. Otherwise,
        all ids are considered. Any ids in skip_id will be silently ignored.

        Filtering on (ra,dec) to determine visibility is currently not very robust and assumes
        that ra_center is close to zero and that subtracting 360 from any ra > 180 is sufficient
        for interval tests.

        Args:
            survey(:class:`descwl.survey.Survey`): Survey parameters used to determine which
                entries are visible.
            render_options(:class:`descwl.render.Options`): Rendering options used to determine
                which entries are visible.

        Yields:
            tuple: Tuple (entry,dx,dy) of the current catalog entry (:class:`astropy.table.Row`)
                and the x,y offsets (:class:`float`) of this entry's centroid from the image
                center in pixels.
        """
        # Calculate the margin size in arcsecs.
        margin_size = 0. if render_options.no_margin else render_options.truncate_radius
        # Calculate the RA,DEC limits of visible entries in degrees.
        arcsec2deg = 1./3600.
        # Calculate scaling between RA and vertical angular separations in the image.
        ra_scale = math.cos(math.radians(self.dec_center))
        ra_size = (0.5*survey.image_width*survey.pixel_scale + margin_size)*arcsec2deg/ra_scale
        ra_min = self.ra_center - ra_size
        ra_max = self.ra_center + ra_size
        dec_size = (0.5*survey.image_height*survey.pixel_scale + margin_size)*arcsec2deg
        dec_min = self.dec_center - dec_size
        dec_max = self.dec_center + dec_size
        # Iterate over all catalog entries.
        for entry in self.table:
            if self.only_id and entry['id'] not in self.only_id:
                continue
            if self.skip_id and entry['id'] in self.skip_id:
                continue
            ra,dec = entry['ra'],entry['dec']
            if ra > 180:
                ra -= 360.
            if ra < ra_min or ra > ra_max or dec < dec_min or dec > dec_max:
                continue
            # If we get this far, the entry is visible.
            dx_arcsecs = (ra - self.ra_center)/arcsec2deg*ra_scale
            dy_arcsecs = (dec - self.dec_center)/arcsec2deg
            yield entry,dx_arcsecs,dy_arcsecs
        
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
        parser.add_argument('--star_catalog-name', type = str, default = None, metavar = 'NAME',
            help = 'Name of catalog file, which must exist and be readable.')
        
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
                