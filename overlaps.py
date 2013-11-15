#!/usr/bin/env python
#######################################################################################
## Created by David Kirkby, University of California, Irvine <dkirkby@uci.edu>
#######################################################################################

import sys
import os
import math
import argparse
import operator

import logging
import galsim
import numpy

def main():

    # Parse command-line args
    parser = argparse.ArgumentParser()
    parser.add_argument("-v", "--verbose", action = "store_true",
        help = "provide more verbose output")
    parser.add_argument("-i", "--input", default = 'catout',
        help = "galsimcat base name of files to read")
    args = parser.parse_args()

    # Configure the GalSim logger
    logging.basicConfig(format="%(message)s", level=logging.INFO, stream=sys.stdout)
    logger = logging.getLogger("overlaps")
    logger.info('Using input prefix %r' % args.input)

    # Scan the group list for groups with more than one member
    overlapGroups = [ ]
    for line in open(args.input + '_groups.dat','r'):
        fields = line.split()
        groupID = int(fields[0])
        groupSize = int(fields[1])
        if groupSize > 1:
            overlapGroups.append(groupID)

    # Get info on the members of each overlap group
    overlaps = { }
    for line in open(args.input + '_catalog.dat','r'):
        fields = line.split()
        flux = float(fields[4])
        z = float(fields[9])
        groupID = int(fields[11])
        if groupID in overlapGroups:
            overlap = overlaps.get(groupID,[ ])
            overlap.append((flux,z))
            overlaps[groupID] = overlap
    print 'Found %d overlap groups' % len(overlaps)

    # Loop over overlap groups
    dzmax = 4.05
    dzmin = -4.05
    binsize = 0.1
    nbins = int(math.ceil((dzmax-dzmin)/binsize))
    histo = numpy.zeros(nbins)
    for (groupID,overlap) in overlaps.iteritems():
        # Get the total flux in this overlap
        sumflux = sum([flux for (flux,z) in overlap])
        # Find the member with the most flux
        (maxflux,maxz) = max(overlap, key = operator.itemgetter(0))
        for (flux,z) in overlap:
            dz = z - maxz
            bin = int(math.floor((dz - dzmin)/binsize))
            if bin >= 0 and bin < nbins:
                histo[bin] += flux/sumflux

    out = open('hist.dat','w')
    for bin in range(nbins):
        dz = dzmin + (bin+0.5)*binsize
        print >>out, '%.2f %f' % (dz,histo[bin]/len(overlaps))
    out.close()

if __name__ == "__main__":
    main()
