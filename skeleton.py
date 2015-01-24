#!/usr/bin/env python
"""Skeleton program to demonstrate reading and analyzing simulation output.

This program reads a simulation output file 'demo.fits' and then loops over
all overlapping groups with exactly two members, with some additional cuts on
the galaxy properties, finally saving images of each pair to an output
file 'pairs.fits'.
"""

import numpy as np
import galsim
import descwl

# Initialize an array of GalSim images that we will save.
images_to_save = [ ]

# Load the results of a previous simulation.
results = descwl.output.Reader('demo.fits').results

# Select the brightest sources from groups with exactly 2 members.
selected = results.select('grp_size==2','grp_rank==0')

# Loop over these pair groups.
for index in selected:
	info = results.table[index]
	# Skip groups where the overlap has essentially no effect on the brighter galaxy.
	if info['purity'] > 0.95: continue
	# Select all (both) members of this group.
	group = results.select('grp_id==%d' % info['db_id'])
	# Create and save an image of just this group.
	image = results.get_subimage(group)
	images_to_save.append(image)

print 'Saving %d images.' % len(images_to_save)
galsim.fits.writeMulti(images_to_save,file_name = 'pairs.fits')
