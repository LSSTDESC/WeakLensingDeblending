#!/usr/bin/env python
"""Fast image simulation using GalSim for weak lensing studies.
"""

import argparse
import os
import os.path

import astropy.io.fits

import descwl

def main():
    # Initialize and parse command-line arguments.
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--verbose', action = 'store_true',
        help = 'Provide verbose output.')
    parser.add_argument('--survey-defaults', action = 'store_true',
        help = 'Print survey camera and observing parameter defaults and exit.')
    catalog_group = parser.add_argument_group('Catalog input',
        'Specify an input catalog of source parameters for simulation.')
    descwl.catalog.Reader.add_args(catalog_group)
    survey_group = parser.add_argument_group('Survey parameters',
        'Specify survey camera and observing parameters.')
    descwl.survey.Survey.add_args(survey_group)
    model_group = parser.add_argument_group('Source model options',
        'Specify options for building source models from catalog parameters.')
    descwl.model.GalaxyBuilder.add_args(model_group)
    render_group = parser.add_argument_group('Model rendering options',
        'Specify options for rendering models as simulated survey observations.')
    descwl.render.Engine.add_args(render_group)
    output_group = parser.add_argument_group('Output control',
        'Specify options to control simulation output.')
    output_group.add_argument('--output-name', default = None, metavar = 'FILE',
        help = 'Base name of FITS output files to write. The ".fits" extension can be omitted.')
    output_group.add_argument('--output-no-clobber', action = 'store_true',
        help = 'Do not overwrite any existing output file with the same name.')
    args = parser.parse_args()

    if args.survey_defaults:
        descwl.survey.Survey.print_defaults()
        return 0

    try:

        catalog = descwl.catalog.Reader.from_args(args)
        if args.verbose:
            print 'Read %d catalog entries from %s' % (len(catalog.table),catalog.catalog_name)

        survey = descwl.survey.Survey.from_args(args)
        if args.verbose:
            print survey.description()

        render_engine = descwl.render.Engine.from_args(survey,args)
        galaxy_builder = descwl.model.GalaxyBuilder.from_args(survey,args)

        hdu_list = None
        if args.output_name:
            name,extension = os.path.splitext(args.output_name)
            if not extension:
                args.output_name += '.fits'
            elif extension.lower() != '.fits':
                print 'Got unexpected output-name extension "%s".' % extension
                return -1
            if not os.access(args.output_name,os.W_OK):
                print 'Requested output file is not writeable: %s' % args.output_name
                return -1
            if args.verbose:
                print 'Simulation output will be saved to %s' % args.output_name
            primary_hdu = astropy.io.fits.PrimaryHDU(data = survey.image.array)
            hdu_list = astropy.io.fits.HDUList([primary_hdu])

        for entry,dx,dy in catalog.visible_entries(survey,render_engine):

            galaxy = galaxy_builder.from_catalog(entry,dx,dy,survey.filter_band)
            if galaxy is None:
                continue

            galaxy_stamps = render_engine.render_galaxy(galaxy)
            if galaxy_stamps is None:
                continue
            if hdu_list:
                data_cube = astropy.io.fits.ImageHDU(data = galaxy_stamps)
                hdu_list.append(data_cube)

        if hdu_list:
            hdu_list.writeto(args.output_name,clobber = not args.output_no_clobber)

    except RuntimeError,e:
        print str(e)
        return -1

if __name__ == '__main__':
    main()
