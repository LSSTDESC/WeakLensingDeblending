#!/usr/bin/env python
"""Fast image simulation using GalSim for weak lensing studies.
"""

import argparse

import descwl

def main():
    # Initialize and parse command-line arguments.
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--verbose', action = 'store_true',
        help = 'Provide verbose output.')
    parser.add_argument('--survey-name', choices = ['LSST','DES','CFHT'], default = 'LSST',
        help = 'Use default camera and observing parameters appropriate for this survey.')
    parser.add_argument('--filter-band', choices = ['u','g','r','i','z','y'], default = 'i',
        help = 'LSST imaging band to simulate')
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
            print 'Simulating %s %s-band survey with %r' % (
                args.survey_name,args.filter_band,survey.args)

        render_engine = descwl.render.Engine.from_args(args)
        galaxy_builder = descwl.model.GalaxyBuilder.from_args(survey,args)

        for entry,dx,dy in catalog.visible_entries(survey,render_engine):

            galaxy_model = galaxy_builder.from_catalog(entry,dx,dy,args.filter_band)
            if galaxy_model is None:
                continue

            galaxy_stamps = render_engine.render_galaxy(galaxy_model,survey)
            if galaxy_stamps is None:
                continue

    except RuntimeError,e:
        print str(e)
        return -1

if __name__ == '__main__':
    main()
