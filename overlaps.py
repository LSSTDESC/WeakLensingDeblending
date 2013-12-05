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
    print 'Found %d overlap groups' % len(overlapGroups)

    # Get info on the members of each overlap group
    overlaps = { }
    ncatalog = 0
    for line in open(args.input + '_catalog.dat','r'):
        ncatalog += 1
        fields = line.split()
        flux = float(fields[4])
        z = float(fields[14])
        groupID = int(fields[16])
        if groupID in overlapGroups:
            overlap = overlaps.get(groupID,[ ])
            overlap.append((flux,z))
            overlaps[groupID] = overlap
    print 'Found %d overlap groups from %d galaxies' % (len(overlaps),ncatalog)

    # Load the full field image and individual stamps
    field = galsim.fits.read(args.input + '_field.fits')
    for i in range(ncatalog):
        stamps = galsim.fits.readCube(args.input + '_stamps.fits',hdu=i+1)
        print i,len(stamps),stamps[0].array.shape

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
