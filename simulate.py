#!/usr/bin/env python
"""Fast image simulation using GalSim for weak lensing studies.
"""

import argparse

import descwl

def main():
    # Initialize and parse command-line arguments.
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    catalog_group = parser.add_argument_group('Catalog input',
        'Specify an input catalog of source parameters for simulation.')
    descwl.catalog.Reader.add_args(catalog_group)
    args = parser.parse_args()

    catalog = descwl.catalog.Reader.from_args(args)
    entries = 0
    for entry in catalog:
        if entries < 10:
            print entry
        entries += 1

if __name__ == '__main__':
    main()
