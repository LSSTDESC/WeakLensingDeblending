#!/usr/bin/env python

"""Create plots to illustrate galaxy parameter error estimation using Fisher matrices.
"""

import argparse

import matplotlib.pyplot as plt

import descwl

def main():
    # Initialize and parse command-line arguments.
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--verbose', action = 'store_true',
        help = 'Provide verbose output.')
    descwl.output.Reader.add_args(parser)
    parser.add_argument('--no-display', action = 'store_true',
        help = 'Do not display the image on screen.')
    parser.add_argument('-o','--output-name',type = str, default = None, metavar = 'FILE',
        help = 'Name of the output file to write.')

    args = parser.parse_args()
    if args.no_display and not args.output_name:
        print 'No display our output requested.'
        return 0

    # Load the analysis results file we will get partial derivative images from.
    try:
        reader = descwl.output.Reader.from_args(args)
        results = reader.results
        if args.verbose:
            print results.survey.description()
    except RuntimeError,e:
        print str(e)
        return -1

    figure = plt.figure()

    if args.output_name:
        figure.savefig(args.output_name)

    if not args.no_display:
        plt.show()

if __name__ == '__main__':
    main()
