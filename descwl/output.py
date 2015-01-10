"""Configure and handle simulation output.

There is a separate :doc:`output page </output>` with details on what goes into the
output and how it is formatted.
"""

import os
import os.path
import inspect

import numpy as np

import astropy.io.fits
import astropy.table

import galsim

import descwl.survey
import descwl.analysis

class Reader(object):
    """Simulation output reader.

    The reader loads the files contents into memory in the constructor and makes them
    available via a `results` data member of type :class:`descwl.analysis.OverlapResults`.

    Args:
        input_name(str): Base name of FITS output files to write. The ".fits" extension
            can be omitted.

    Raises:
        RuntimeError: Unable to initialize FITS input file.
    """
    def __init__(self,input_name):
        if not input_name:
            raise RuntimeError('Missing required input-name parameter.')
        self.input_name = input_name
        name,extension = os.path.splitext(self.input_name)
        if not extension:
            self.input_name += '.fits'
        elif extension.lower() != '.fits':
            raise RuntimeError('Got unexpected input-name extension "%s".' % extension)
        # Without memmap=False, reading all the postage stamp HDUs crashes with
        # "OSError: [Errno 24] Too many open files" since each access to hdu_list[n].data
        # adds a new open file descriptor.
        try:
            self.hdu_list = astropy.io.fits.open(self.input_name,mode='readonly',memmap=False)
        except IOError,e:
            raise RuntimeError(str(e))
        # Reconstruct the survey object for these results.
        header = self.hdu_list[0].header
        survey_args = { }
        for parameter_name in descwl.survey.Survey._parameter_names:
            # Fits header keyword names are truncated at length 8.
            survey_args[parameter_name] = header[parameter_name[:8].upper()]
        survey = descwl.survey.Survey(**survey_args)
        # Load the simulated image into the survey object.
        image_data = self.hdu_list[0].data
        survey.image.array[:] = image_data
        # Passing an HDUList to Table.read does not seem to be documented but works as expected.
        table = astropy.table.Table.read(self.hdu_list,hdu=1)
        # Load individual stamps and reconstruct the corresponding bounds objects.
        stamps,bounds = [ ],[ ]
        num_galaxies = len(table)
        for index in range(num_galaxies):
            # Stamp HDUs start after the full image HDU[0] and the analysis table HDU[1].
            hdu = self.hdu_list[index+2]
            datacube = np.copy(hdu.data)
            stamps.append(datacube)
            nstamps,height,width = datacube.shape
            x_min,y_min = hdu.header['X_MIN'],hdu.header['Y_MIN']
            bounds.append(galsim.BoundsI(x_min,x_min+width-1,y_min,y_min+height-1))
        # Save file contents as a results object.
        self.results = descwl.analysis.OverlapResults(survey,table,stamps,bounds)
        self.hdu_list.close()

    @staticmethod
    def add_args(parser):
        """Add command-line arguments for constructing a new :class:`Reader`.

        The added arguments are our constructor parameters with '_' replaced by '-' in the names.

        Args:
            parser(argparse.ArgumentParser): Arguments will be added to this parser object using its
                add_argument method.
        """
        parser.add_argument('-i','--input-name', default = None, metavar = 'FILE',
            help = 'Base name of FITS output files to read. The ".fits" extension can be omitted.')

    @classmethod
    def from_args(cls,args):
        """Create a new :class:`Reader` object from a set of arguments.

        Args:
            survey(descwl.survey.Survey): Simulated survey to describe with FITS header keywords.
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

class Writer(object):
    """Simulation output writer.

    See the :doc:`output page </output>` for details on the output contents and formatting.

    Args:
        survey(descwl.survey.Survey): Simulated survey to describe with FITS header keywords.
        output_name(str): Base name of FITS output files to write. The ".fits" extension
            can be omitted.
        no_catalog(bool): Do not save the results catalog.
        no_stamps(bool): Do not save postage stamp datacubes for each source.
        output_no_clobber(bool): Do not overwrite any existing output file with the same name.

    Raises:
        RuntimeError: Unable to initialize FITS output file.
    """
    def __init__(self,survey,output_name,no_stamps,no_catalog,output_no_clobber):
        self.survey = survey
        self.output_name = output_name
        self.no_stamps = no_stamps
        self.no_catalog = no_catalog
        self.output_no_clobber = output_no_clobber
        self.hdu_list = None
        if self.output_name:
            name,extension = os.path.splitext(self.output_name)
            if not extension:
                self.output_name += '.fits'
            elif extension.lower() != '.fits':
                raise RuntimeError('Got unexpected output-name extension "%s".' % extension)
            try:
                with open(self.output_name,'rb') as file:
                    # File already exists. Can we clobber it?
                    if self.output_no_clobber:
                        raise RuntimeError('Cannot clobber existing file %r' %
                            self.output_name)
            except IOError:
                pass
            # Try to open an empty file, deleting any previous contents.
            try:
                with open(self.output_name,'wb') as file:
                    pass
            except IOError:
                raise RuntimeError('File is not writeable %r' % args.output_name)
            # Open the now-empty file as a FITS file.
            self.hdu_list = astropy.io.fits.open(self.output_name, mode = 'update',
                memmap = False)
            assert len(self.hdu_list) == 0, 'Expected new FITS file to be empty.'
            # Add an empty primary HDU for now.
            primary_hdu = astropy.io.fits.PrimaryHDU()
            self.hdu_list.append(primary_hdu)

    def description(self):
        """Describe our output configuration.

        Returns:
            str: Description of the the rendering configuration.
        """
        return 'Simulation output will be saved to %s' % self.output_name

    def finalize(self,results):
        """Save analysis results and close the output file, if any.

        Args:
            :class:`descwl.analysis.OverlapResults`: Overlap analysis results.
        """
        if self.hdu_list is None:
            return
        # Fill in the primary HDU image data and headers.
        primary_hdu = self.hdu_list[0]
        primary_hdu.data = results.survey.image.array
        for key,value in results.survey.args.iteritems():
            primary_hdu.header.set(key[:8],value)
        if not self.no_catalog:
            # Save the analysis results table in HDU[1].
            table = astropy.io.fits.BinTableHDU.from_columns(np.array(results.table))
            self.hdu_list.append(table)
        if not self.no_stamps:
            # Save each stamp datacube.
            for stamps,bounds in zip(results.stamps,results.bounds):
                data_cube = astropy.io.fits.ImageHDU(data = stamps)
                data_cube.header['X_MIN'] = bounds.xmin
                data_cube.header['Y_MIN'] = bounds.ymin
                self.hdu_list.append(data_cube)
        # Write and close our FITS file.
        self.hdu_list.close()

    @staticmethod
    def add_args(parser):
        """Add command-line arguments for constructing a new :class:`Writer`.

        The added arguments are our constructor parameters with '_' replaced by '-' in the names.

        Args:
            parser(argparse.ArgumentParser): Arguments will be added to this parser object using its
                add_argument method.
        """
        parser.add_argument('--output-name', default = None, metavar = 'FILE',
            help = 'Base name of FITS output files to write. The ".fits" extension can be omitted.')
        parser.add_argument('--no-catalog', action = 'store_true',
            help = 'Do not save the results catalog.')
        parser.add_argument('--no-stamps', action = 'store_true',
            help = 'Do not save postage stamp datacubes for each source.')
        parser.add_argument('--output-no-clobber', action = 'store_true',
            help = 'Do not overwrite any existing output file with the same name.')

    @classmethod
    def from_args(cls,survey,args):
        """Create a new :class:`Writer` object from a set of arguments.

        Args:
            survey(descwl.survey.Survey): Simulated survey to describe with FITS header keywords.
            args(object): A set of arguments accessed as a :py:class:`dict` using the
                built-in :py:func:`vars` function. Any extra arguments beyond those defined
                in :func:`add_args` will be silently ignored.

        Returns:
            :class:`Writer`: A newly constructed Reader object.
        """
        # Look up the named constructor parameters.
        pnames = (inspect.getargspec(cls.__init__)).args[1:]
        # Get a dictionary of the arguments provided.
        args_dict = vars(args)
        # Filter the dictionary to only include constructor parameters.
        filtered_dict = { key:args_dict[key] for key in (set(pnames) & set(args_dict)) }
        return cls(survey,**filtered_dict)
