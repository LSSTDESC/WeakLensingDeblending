"""Configure and handle simulation output.

There is a separate :doc:`output page </output>` with details on what goes into the
output and how it is formatted.
"""
from __future__ import print_function, division

import os
import os.path
import inspect

import numpy as np

import astropy.table
import astropy.io.fits

import fitsio

import galsim

import descwl.survey
import descwl.analysis

class Reader(object):
    """Simulation output reader.

    The reader loads the file's contents into memory in the constructor and makes them available
    via a `results` data member of type :class:`descwl.analysis.OverlapResults`. Loading of stamp
    datacubes can be defered until stamps are acutally accessed using the defer_stamp_loading
    argument below.

    We use `fitsio <https://github.com/esheldon/fitsio>`_ to read files since it performs
    significantly faster than :mod:`astropy.io.fits` for files with many (~100K) relatively small
    HDUs (although :mod:`astropy.io.fits` outperforms `fitsio` for writing these files).

    Args:
        input_name(str): Base name of FITS output files to write. The ".fits" extension
            can be omitted.
        defer_stamp_loading(bool): Defer the loading of stamp datacubes until they are actually
            accessed via :func:`descwl.OverlapResults.get_stamp` of our `results` data
            member. This speeds up the constructor substantially and reduces memory
            requirements when relatively few stamps are actually needed.

    Raises:
        RuntimeError: Unable to initialize FITS input file.
    """
    def __init__(self,input_name,defer_stamp_loading = True):
        if not input_name:
            raise RuntimeError('Missing required input-name parameter.')
        self.input_name = input_name
        name,extension = os.path.splitext(self.input_name)
        if not extension:
            self.input_name += '.fits'
        elif extension.lower() != '.fits':
            raise RuntimeError('Got unexpected input-name extension "%s".' % extension)
        try:
            self.fits = fitsio.FITS(self.input_name,mode = fitsio.READONLY)
        except ValueError as e:
            raise RuntimeError(str(e))
        # Reconstruct the survey object for these results.
        header = self.fits[0].read_header()
        num_slices = header['NSLICES']
        survey_args = { }
        for parameter_name in descwl.survey.Survey._parameter_names:
            # Fits header keyword names are truncated at length 8.
            value = header[parameter_name[-8:].upper()]
            # String values are padded on the right with spaces.
            survey_args[parameter_name] = value.rstrip() if type(value) is str else value
        survey = descwl.survey.Survey(**survey_args)
        survey.psf_sigma_m = header['PSF_SIGM']
        survey.psf_sigma_p = header['PSF_SIGP']
        survey.psf_size_hsm = header['PSF_HSM']
        # Load the simulated image into the survey object.
        survey.image.array[:] = self.fits[0].read()

        table = None
        stamp_hdu_offset = 1
        if len(self.fits) > 1 and type(self.fits[1]) is fitsio.fitslib.TableHDU:
            table = astropy.table.Table(self.fits[1].read(),copy = False)
            stamp_hdu_offset += 1

        stamps,bounds = [ ],[ ]
        if len(self.fits) > stamp_hdu_offset:
            if table is None:
                raise RuntimeError('Missing required table for reconstructing bounding boxes.')
            # Load individual stamps and reconstruct the corresponding bounds objects.
            for hdu_index in range(stamp_hdu_offset,len(self.fits)):
                # Reconstruct the stamp bounds from the corresponding row of the table.
                row = table[hdu_index - stamp_hdu_offset]
                bounds.append(galsim.BoundsI(row['xmin'],row['xmax'],row['ymin'],row['ymax']))
                # Load the stamp datacube, now or later.
                if defer_stamp_loading:
                    # Make sure we bind the current value of hdu_index, not the variable itself.
                    stamps.append(lambda index=hdu_index: self.fits[index].read())
                else:
                    stamps.append(self.fits[hdu_index].read())
        # Save file contents as a results object.
        self.results = descwl.analysis.OverlapResults(survey,table,stamps,bounds,num_slices)
        if not defer_stamp_loading:
            self.fits.close()

    @staticmethod
    def add_args(parser):
        """Add command-line arguments for constructing a new :class:`Reader`.

        The added arguments are our constructor arguments with '_' replaced by '-' in the names.
        The defer_stamp_loading constructor argument is not added here.

        Args:
            parser(argparse.ArgumentParser): Arguments will be added to this parser object using its
                add_argument method.
        """
        parser.add_argument('-i','--input-name', default = None, metavar = 'FILE',
            help = 'Base name of FITS output files to read. The ".fits" extension can be omitted.')

    @classmethod
    def from_args(cls,defer_stamp_loading,args):
        """Create a new :class:`Reader` object from a set of arguments.

        Args:
            defer_stamp_loading(bool): Defer the loading of stamp datacubes until they are actually
                accessed via :func:`descwl.OverlapResults.get_stamp` of our `results` data
                member. This speeds up the constructor substantially and reduces memory
                requirements when relatively few stamps are actually needed.
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
        return cls(defer_stamp_loading = defer_stamp_loading,**filtered_dict)

class Writer(object):
    """Simulation output writer.

    See the :doc:`output page </output>` for details on the output contents and formatting.

    We use :mod:`astropy.io.fits` to write files since it performs significanly faster than
    `fitsio <https://github.com/esheldon/fitsio>`_ for files with many (~100K) relatively small
    HDUs (although `fitsio` outperforms :mod:`astropy.io.fits` for reading these files).
    See `this issue <https://github.com/esheldon/fitsio/issues/32>`_ for details and updates.

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
        self.fits = None
        if self.output_name:
            name,extension = os.path.splitext(self.output_name)
            if not extension:
                self.output_name += '.fits'
            elif extension.lower() != '.fits':
                raise RuntimeError('Got unexpected output-name extension "%s".' % extension)
            if not self.output_no_clobber:
                # Remove the file if it already exists.
                try:
                    os.unlink(self.output_name)
                except OSError:
                    pass
            # Try to open a new file.
            try:
                self.fits = astropy.io.fits.open(self.output_name,mode = 'ostream',memmap = False)
                self.fits.append(astropy.io.fits.PrimaryHDU())
            except IOError as e:
                raise RuntimeError(str(e))

    def description(self):
        """Describe our output configuration.

        Returns:
            str: Description of the the rendering configuration.
        """
        return 'Simulation output will be saved to %s' % self.output_name

    def finalize(self,results,trace):
        """Save analysis results and close the output file, if any.

        This call builds the entire FITS file in memory then writes it to disk.

        Args:
            :class:`descwl.analysis.OverlapResults`: Overlap analysis results.
            trace(callable): Function to call for tracing resource usage. Will be
                called with a brief :class:`str` description of each checkpoint.
        """
        trace('Writer.finalize begin')
        if self.fits is None:
            return
        # Write the simulated survey image into the primary HDU.
        hdu = self.fits[0]
        hdu.data = results.survey.image.array
        # Copy our Survey ctor args into the primary HDU header.
        hdu.header['NSLICES'] = results.num_slices
        hdu.header['PSF_SIGM'] = self.survey.psf_sigma_m
        hdu.header['PSF_SIGP'] = self.survey.psf_sigma_p
        hdu.header['PSF_HSM'] = self.survey.psf_size_hsm
        for key,value in results.survey.args.iteritems():
            # Fits keyword headers are truncated at length 8. We use the last 8 chararacters
            # to ensure that they are unique.
            hdu.header[key[-8:]] = value
        trace('wrote primary hdu')
        if not self.no_catalog:
            # Save the analysis results table in HDU[1].
            hdu = astropy.io.fits.BinTableHDU.from_columns(np.array(results.table))
            self.fits.append(hdu)
            trace('wrote table')
        if not self.no_stamps:
            # Save each stamp datacube.
            for cube in results.stamps:
                hdu = astropy.io.fits.ImageHDU(data = cube)
                self.fits.append(hdu)
                trace('wrote datacube')
        # Write and close our FITS file. Nothing is actually written to disk until we call flush().
        self.fits.flush()
        self.fits.close()
        trace('Writer.finalize end')

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
