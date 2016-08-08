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
catalog = results.table

# Select the brightest sources from groups with exactly 2 members.
selected = results.select('grp_size==2','grp_rank==0')

# Loop over these pair groups.
for index in selected:
	# Skip groups where the overlap has essentially no effect on the brighter galaxy.
	if catalog['purity'][index] > 0.95: continue
	# Select all (both) members of this group.
	grp_id = catalog['db_id'][index]
	group = results.select('grp_id==%ld' % grp_id)
	# Print the fluxes of each galaxy (in order of increasing SNR).
	print 'fluxes for group id %ld are %.3f,%.3f electrons.' % (
		grp_id,catalog['flux'][group[0]],catalog['flux'][group[1]])
	# Create and save an image of just this group.
	image = results.get_subimage(group)
	images_to_save.append(image)
	# Save images of each galaxy individually, using the same bounding box.
	for galaxy_index in group:
		galaxy = image.copy()
		galaxy.array[:] = 0.
		stamp = results.get_stamp(galaxy_index)
		galaxy[stamp.bounds] = stamp
		images_to_save.append(galaxy)

print 'Saving %d images.' % len(images_to_save)
galsim.fits.writeMulti(images_to_save,file_name = 'pairs.fits')
