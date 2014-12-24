#!/usr/bin/env python
"""Display simulated images and analysis results generated by the simulate program.
"""

import argparse
import os.path

import numpy as np
import numpy.ma

import matplotlib.pyplot as plt

import astropy.io.fits

import descwl

def main():
    # Initialize and parse command-line arguments.
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--verbose', action = 'store_true',
        help = 'Provide verbose output.')
    descwl.output.Reader.add_args(parser)
    parser.add_argument('-o','--output-name',type = str, default = None, metavar = 'FILE',
        help = 'Name of the output file to write.')

    select_group = parser.add_argument_group('Object selection options')
    select_group.add_argument('--galaxy', type = int, default = None, metavar = 'ID',
        help = 'Highlight the galaxy with this unique identifier.')
    select_group.add_argument('--crop', action = 'store_true',
        help = 'Crop the displayed pixels around the selected objects.')

    display_group = parser.add_argument_group('Display options')
    display_group.add_argument('--dpi', type = float, default = 64.,
        help = 'Number of pixels per inch to use for display.')
    display_group.add_argument('--magnification', type = float,
        default = 1, metavar = 'MAG',
        help = 'Magnification factor to use for display.')
    display_group.add_argument('--colormap', type = str,
        default = 'gray_r', metavar = 'CMAP',
        help = 'Matplotlib colormap name to use for background pixel values.')
    display_group.add_argument('--highlight', type = str,
        default = 'hot_r', metavar = 'CMAP',
        help = 'Matplotlib colormap name to use for highlighted pixel values.')
    display_group.add_argument('--clip-lo-percentile', type = float,
        default = 0.0, metavar = 'PCT',
        help = 'Clip pixels with values below this percentile for the image.')
    display_group.add_argument('--clip-hi-percentile', type = float,
        default = 90.0, metavar = 'PCT',
        help = 'Clip pixels with values above this percentile for the image.')

    args = parser.parse_args()

    # Load the analysis results file we will display from.
    reader = descwl.output.Reader.from_args(args)
    results = reader.results
    if args.verbose:
        print results.survey.description()

    # Perform object selection
    selected_image = None
    selected_indices = [ ]
    if args.galaxy:
        index = results.find_galaxy(args.galaxy)
        selected_indices.append(index)
        selected_image = results.get_galaxy_image(index)

    # Use the selected pixels to determine the z scaling.
    selected_pixels = selected_image.array
    non_zero_pixels = (selected_pixels > 0)
    vmin,vmax = np.percentile(selected_pixels[non_zero_pixels],
        q = (args.clip_lo_percentile,args.clip_hi_percentile))

    def znorm(pixels):
        return (np.clip(pixels,vmin,vmax) - vmin)/(vmax-vmin)

    # See http://ds9.si.edu/ref/how.html#Scales
    def zscale(pixels):
        return np.sqrt(znorm(pixels))

    # Calculate our viewing bounds.
    if args.crop and selected_image is not None:
        view_bounds = selected_image.bounds
    else:
        view_bounds = results.survey.image.bounds

    # Initialize a matplotlib figure to display our view bounds.
    view_width = view_bounds.xmax - view_bounds.xmin + 1
    view_height = view_bounds.ymax - view_bounds.ymin + 1
    fig_height = args.magnification*(view_height/args.dpi)
    fig_width = args.magnification*(view_width/args.dpi)
    figure = plt.figure(figsize = (fig_width,fig_height),frameon = False,dpi = args.dpi)
    axes = plt.Axes(figure, [0., 0., 1., 1.])
    axes.axis(xmin=view_bounds.xmin,xmax=view_bounds.xmax+1,
        ymin=view_bounds.ymin,ymax=view_bounds.ymax+1)
    axes.set_axis_off()
    figure.add_axes(axes)

    def show_image(image,masked,**kwargs):
        overlap = image.bounds & view_bounds
        xlo = overlap.xmin
        xhi = overlap.xmax + 1
        ylo = overlap.ymin
        yhi = overlap.ymax + 1
        overlap_pixels = image[overlap].array
        z = zscale(overlap_pixels)
        if masked:
            # Only show non-zero pixels.
            z = numpy.ma.masked_where(overlap_pixels == 0,z)
        axes.imshow(z,extent = (xlo,xhi,ylo,yhi),
            aspect = 'equal',origin = 'lower',interpolation = 'nearest',**kwargs)

    # Plot the full simulated image using the background colormap.
    show_image(results.survey.image,masked = False,cmap = args.colormap)

    # Overplot the selected objects showing only non-zero pixels.
    show_image(selected_image,masked = True,cmap = args.highlight)

    # Draw a crosshair at the centroid of selected objects.
    for index in selected_indices:
        # Calculate the centroid position relative to the lower-left corner in pixel units.
        x_pixels = (0.5*results.survey.image_width +
            results.table['dx'][index]/results.survey.pixel_scale)
        y_pixels = (0.5*results.survey.image_height +
            results.table['dy'][index]/results.survey.pixel_scale)
        axes.plot(x_pixels,y_pixels,'w+',markeredgewidth = 2,markersize = 24)

    if args.output_name:
        figure.savefig(args.output_name,dpi = args.dpi)

    plt.show()

if __name__ == '__main__':
    main()
