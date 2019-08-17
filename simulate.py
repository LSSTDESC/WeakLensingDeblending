#!/usr/bin/env python
"""Fast image simulation using GalSim for weak lensing studies.
"""
from __future__ import print_function, division

import argparse

import descwl

def main():
    # Initialize and parse command-line arguments.
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--verbose', action = 'store_true',
        help = 'Provide verbose output.')
    parser.add_argument('--no-analysis', action = 'store_true',
        help = 'Don\'t run analysis.')
    parser.add_argument('--survey-defaults', action = 'store_true',
        help = 'Print survey camera and observing parameter defaults and exit.')
    parser.add_argument('--memory-trace', action = 'store_true',
        help = 'Trace memory usage (requires the psutil module).')
    analysis_group = parser.add_argument_group('Analysis options',
        'Specify analysis options')
    descwl.analysis.OverlapAnalyzer.add_args(analysis_group)
    catalog_group = parser.add_argument_group('Catalog input',
        'Specify an input catalog of source parameters for simulation.')
    descwl.catalog.Reader.add_args(catalog_group)
    catalog_star_group = parser.add_argument_group('Star catalog input',
        'Specify an input star catalog for simulation')
    descwl.catalog.ReaderStar.add_args(catalog_star_group)
    survey_group = parser.add_argument_group('Survey parameters',
        'Specify survey camera and observing parameters.')
    descwl.survey.Survey.add_args(survey_group)
    model_group = parser.add_argument_group('Source model options',
        'Specify options for building source models from catalog parameters.')
    descwl.model.GalaxyBuilder.add_args(model_group)
    model_star_group = parser.add_argument_group('Source model options',
        'Specify options for building source models from catalog parameters.')
    descwl.model.StarBuilder.add_args(model_star_group)
    render_group = parser.add_argument_group('Model rendering options',
        'Specify options for rendering models as simulated survey observations.')
    descwl.render.Engine.add_args(render_group)
    output_group = parser.add_argument_group('Output control',
        'Specify options to control simulation output.')
    descwl.output.Writer.add_args(output_group)
    args = parser.parse_args()

    if args.survey_defaults:
        descwl.survey.Survey.print_defaults()
        return 0

    if args.catalog_name is None and args.star_catalog_name is None: 
        raise RuntimeError("At least one of catalog_name or start_catalog_name must be specifed in order to simulate.")

    if args.no_fisher and args.add_lmfit: 
        raise RuntimeError("Fisher calculation is necessary to run fits with lmfit.")

    if args.no_fisher and args.calculate_bias: 
        raise RuntimeError("Bias calculation requires fisher analysis.")


    trace = descwl.trace.Memory(args.memory_trace)
    trace('begin')

    try:
        if args.catalog_name!=None:
            catalog = descwl.catalog.Reader.from_args(args)
        if args.star_catalog_name!=None:
            star_catalog = descwl.catalog.ReaderStar.from_args(args)
        if args.verbose:
            if args.catalog_name!=None:
                print('Read %d catalog entries from %s' % (len(catalog.table),catalog.catalog_name))
            if args.star_catalog_name!=None:
                print('Read %d catalog entries from %s' % (len(star_catalog.table),star_catalog.star_catalog_name))
        survey = descwl.survey.Survey.from_args(args)
        if args.verbose:
            print(survey.description())
        if args.catalog_name!=None:
            galaxy_builder = descwl.model.GalaxyBuilder.from_args(survey,args)
        if args.star_catalog_name!=None:
            star_builder = descwl.model.StarBuilder.from_args(survey,args)

        render_engine = descwl.render.Engine.from_args(survey,args)
        if args.verbose:
            print(render_engine.description())

        analyzer = descwl.analysis.OverlapAnalyzer(survey,args.no_hsm, not args.add_lmfit, args.no_fisher, args.calculate_bias, args.no_analysis, args.add_noise)

        output = descwl.output.Writer.from_args(survey,args)
        if args.verbose:
            print(output.description())


        trace('initialized')
        if args.catalog_name!=None:
            for entry,dx,dy in catalog.potentially_visible_entries(survey,render_engine):

                try:
                    galaxy = galaxy_builder.from_catalog(entry,dx,dy,survey.filter_band)
                    stamps, bounds = render_engine.render_galaxy(
                        galaxy, args.variations_x, args.variations_s, args.variations_g, args.no_fisher, args.calculate_bias, args.no_analysis)
                    analyzer.add_galaxy(galaxy,stamps,bounds)
                    trace('render')

                except (descwl.model.SourceNotVisible,descwl.render.SourceNotVisible):
                    pass

        if args.star_catalog_name!=None:
            for entry,dx,dy in star_catalog.potentially_visible_entries(survey,render_engine):

                try:
                    star = star_builder.from_catalog(entry,dx,dy,survey.filter_band)
                    stamps,bounds = render_engine.render_star(star, args.variations_x, args.variations_s, args.variations_g, args.no_fisher)
                    analyzer.add_star(star,stamps,bounds)
                    trace('render')

                except (descwl.model.SourceNotVisible,descwl.render.SourceNotVisible):
                    pass
        results = analyzer.finalize(args.verbose,trace)
        output.finalize(results,trace)

    except RuntimeError as e:
        print(str(e))
        return -1

if __name__ == '__main__':
    main()
