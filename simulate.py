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
    descwl.output.Writer.add_args(output_group)
    args = parser.parse_args()

    if args.survey_defaults:
        descwl.survey.Survey.print_defaults()
        return 0

    if args.memory_trace:
        try:
            import psutil
            import os
            this_process = psutil.Process(os.getpid())
            global last_usage
            last_usage = this_process.get_memory_info()[0]
            print 'last_usage',last_usage
            def trace(label):
                global last_usage
                usage = this_process.get_memory_info()[0]
                print '%s memory usage: %.3f Mb (%+d bytes)' % (label,
                    usage/float(2**20),usage-last_usage)
                last_usage = usage
        except:
            print 'Unable to initialize memory-trace. Is psutils installed?'
            return -1
    else:
        def trace(label): pass

    trace('begin')

    try:

        catalog = descwl.catalog.Reader.from_args(args)
        if args.verbose:
            print 'Read %d catalog entries from %s' % (len(catalog.table),catalog.catalog_name)

        survey = descwl.survey.Survey.from_args(args)
        if args.verbose:
            print survey.description()

        galaxy_builder = descwl.model.GalaxyBuilder.from_args(survey,args)

        render_engine = descwl.render.Engine.from_args(survey,args)
        if args.verbose:
            print render_engine.description()

        analyzer = descwl.analysis.OverlapAnalyzer(survey)

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

            except (descwl.model.SourceNotVisible,descwl.render.SourceNotVisible):
                pass

        trace('begin analyze')
        results = analyzer.finalize(args.verbose)
        trace('analyzer.finalize')
        output.finalize(results)
        trace('output.finalize')

    except RuntimeError,e:
        print str(e)
        return -1

if __name__ == '__main__':
    main()
