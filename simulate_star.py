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
    parser.add_argument('--survey-defaults', action = 'store_true',
        help = 'Print survey camera and observing parameter defaults and exit.')
    parser.add_argument('--memory-trace', action = 'store_true',
        help = 'Trace memory usage (requires the psutil module).')
    catalog_group = parser.add_argument_group('Catalog input',
        'Specify an input catalog of source parameters for simulation.')
    descwl.catalog_star.Reader.add_args(catalog_group)
    catalog_star_group = parser.add_argument_group('Star catalog input',
        'Specify an input star catalog for simulation')
    descwl.catalog_star.ReaderStar.add_args(catalog_star_group)
    survey_group = parser.add_argument_group('Survey parameters',
        'Specify survey camera and observing parameters.')
    descwl.survey.Survey.add_args(survey_group)
    model_group = parser.add_argument_group('Source model options',
        'Specify options for building source models from catalog parameters.')
    descwl.model_star.GalaxyBuilder.add_args(model_group)
    model_star_group = parser.add_argument_group('Source model options',
        'Specify options for building source models from catalog parameters.')
    descwl.model_star.StarBuilder.add_args(model_star_group)
    render_group = parser.add_argument_group('Model rendering options',
        'Specify options for rendering models as simulated survey observations.')
    descwl.render_star.Engine.add_args(render_group)
    output_group = parser.add_argument_group('Output control',
        'Specify options to control simulation output.')
    descwl.output.Writer.add_args(output_group)
    args = parser.parse_args()

    if args.survey_defaults:
        descwl.survey.Survey.print_defaults()
        return 0

    trace = descwl.trace.Memory(args.memory_trace)
    trace('begin')

    try:

        catalog = descwl.catalog_star.Reader.from_args(args)
        star_catalog = descwl.catalog_star.ReaderStar.from_args(args)
        if args.verbose:
            print 'Read %d catalog entries from %s' % (len(catalog.table),catalog.catalog_name)
            print 'Read %d catalog entries from %s' % (len(star_catalog.table),star_catalog.star_catalog_name)
        survey = descwl.survey.Survey.from_args(args)
        if args.verbose:
            print survey.description()

        galaxy_builder = descwl.model_star.GalaxyBuilder.from_args(survey,args)
		
        star_builder = descwl.model_star.StarBuilder.from_args(survey,args)
        
        render_engine = descwl.render_star.Engine.from_args(survey,args)
        if args.verbose:
            print render_engine.description()

        analyzer = descwl.analysis_star.OverlapAnalyzer(survey)

        output = descwl.output.Writer.from_args(survey,args)
        if args.verbose:
            print output.description()

        trace('initialized')

        for entry,dx,dy in catalog.potentially_visible_entries(survey,render_engine):

            try:
                galaxy = galaxy_builder.from_catalog(entry,dx,dy,survey.filter_band)
                stamps,bounds = render_engine.render_galaxy(galaxy)
                analyzer.add_galaxy(galaxy,stamps,bounds)
                trace('render')

            except (descwl.model_star.SourceNotVisible,descwl.render_star.SourceNotVisible):
                pass

        
        for entry,dx,dy in star_catalog.potentially_visible_entries(survey,render_engine):

            try:
                star = star_builder.from_catalog(entry,dx,dy,survey.filter_band)
                stamps,bounds = render_engine.render_star(star)
                analyzer.add_star(star,stamps,bounds)
                trace('render')

            except (descwl.model_star.SourceNotVisible,descwl.render_star.SourceNotVisible):
                pass

        results = analyzer.finalize(args.verbose,trace)
        output.finalize(results,trace)

    except RuntimeError,e:
        print str(e)
        return -1

if __name__ == '__main__':
    main()
