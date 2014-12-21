"""Configure and handle simulation output.

There is a separate :doc:`output page </output>` with details on what goes into the
output and how it is formatted.
"""

import os
import os.path
import inspect

import numpy as np

import astropy.io.fits

import descwl.survey

class Reader(object):
    """Simulation output reader.

    Args:
        input_name(str): Base name of FITS output files to write. The ".fits" extension
            can be omitted.

    Raises:
        RuntimeError: Unable to initialize FITS input file.
    """
    def __init__(self,input_name):
        self.input_name = input_name
        name,extension = os.path.splitext(self.input_name)
        if not extension:
            self.input_name += '.fits'
        elif extension.lower() != '.fits':
            raise RuntimeError('Got unexpected input-name extension "%s".' % extension)
        self.hdu_list = astropy.io.fits.open(self.input_name)
        # Reconstruct the survey object for these results.
        header = self.hdu_list[0].header
        survey_args = { }
        for parameter_name in descwl.survey.Survey._parameter_names:
            # Fits header keyword names are truncated at length 8.
            survey_args[parameter_name] = header[parameter_name[:8].upper()]
        self.survey = descwl.survey.Survey(**survey_args)
        # Load the simulated image into the survey object.
        image_data = self.hdu_list[0].data
        self.survey.image.array[:] = image_data

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
        output_no_clobber(bool): Do not overwrite any existing output file with the same name.

    Raises:
        RuntimeError: Unable to initialize FITS output file.
    """
    def __init__(self,survey,output_name,output_no_clobber):
        self.output_name = output_name
        self.output_no_clobber = output_no_clobber
        self.hdu_list = None
        if self.output_name:
            name,extension = os.path.splitext(self.output_name)
            if not extension:
                self.output_name += '.fits'
            elif extension.lower() != '.fits':
                raise RuntimeError('Got unexpected output-name extension "%s".' % extension)
            if not os.access(self.output_name,os.W_OK):
                raise RuntimeError('Requested output file is not writeable: %s' % self.output_name)
            primary_hdu = astropy.io.fits.PrimaryHDU(data = survey.image.array)
            for key,value in survey.args.iteritems():
                primary_hdu.header.set(key[:8],value)
            self.hdu_list = astropy.io.fits.HDUList([primary_hdu])

    def description(self):
        """Describe our output configuration.

        Returns:
            str: Description of the the rendering configuration.
        """
        return 'Simulation output will be saved to %s' % self.output_name

    def save_stamps(self,stamps,bounds):
        """Save a datacube of postage stamp images for a single source.

        Args:
            stamps(:class:`numpy.ndarray`): Array of shape (nstamp,width,height) containing
                pixel values for nstamp stamps of dimensions (width,height).
            bounds(galsim.BoundsI): Bounds of the stamps in the full simulated survey image.
        """
        if self.hdu_list:
            data_cube = astropy.io.fits.ImageHDU(data = stamps)
            data_cube.header['X_MIN'] = bounds.xmin
            data_cube.header['Y_MIN'] = bounds.ymin
            self.hdu_list.append(data_cube)

    def finalize(self,results):
        """Save analysis results and close the output file, if any.

        Args:
            :class:`astropy.table.Table`: Table of analysis results with one row per galaxy.
        """
        if self.hdu_list:
            table = astropy.io.fits.BinTableHDU.from_columns(np.array(results))
            self.hdu_list.insert(1,table)
            self.hdu_list.writeto(self.output_name,clobber = not self.output_no_clobber)

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
