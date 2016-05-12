#!/usr/bin/env python

"""Create plots to illustrate galaxy parameter error estimation using Fisher matrices.
"""

import argparse

import numpy as np

import matplotlib.pyplot as plt

import astropy.table

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
    parser.add_argument('--galaxy', type = int,
        default = None, metavar = 'ID',
        help = 'Use the galaxy with this database ID, ignoring any overlaps.')
    parser.add_argument('--group', type = int,
        default = None, metavar = 'ID',
        help = 'Use the overlapping group of galaxies with this group ID.')
    parser.add_argument('--partials', action = 'store_true',
        help = 'Show partial derivative images (instead of Fisher matrix images).')
    parser.add_argument('--matrix', action = 'store_true',
        help = 'Show summed Fisher matrix elements (instead of Fisher matrix images).')
    parser.add_argument('--covariance', action = 'store_true',
        help = 'Show covariance matrix elements (instead of Fisher matrix images).')
    parser.add_argument('--correlation', action = 'store_true',
        help = 'Show correlation matrix elements (instead of Fisher matrix images).')
    parser.add_argument('--bias', action = 'store_true',
        help = 'Show bias of each of the selected galaxies parameters.')

    display_group = parser.add_argument_group('Display options')
    display_group.add_argument('--figure-size', type = float,
        default = 12., metavar = 'SIZE',
        help = 'Size of the longest edge of the figure in inches.')
    display_group.add_argument('--colormap', type = str,
        default = 'RdBu', metavar = 'CMAP',
        help = 'Matplotlib colormap name to use.')
    display_group.add_argument('--no-labels', action = 'store_true',
        help = 'Do not display any text labels.')
    display_group.add_argument('--label-color', type = str,
        default = 'greenyellow', metavar = 'COL',
        help = 'Matplotlib color name to use for label text.')
    display_group.add_argument('--label-size', type = str,
        default = 'medium', metavar = 'SIZE',
        help = 'Matplotlib font size specification in points or relative (small,large,...)')
    display_group.add_argument('--value-format', type = str,
        default = '%.3g', metavar = 'FMT',
        help = 'Printf format to use for matrix element values.')
    display_group.add_argument('--clip-percentile', type = float,
        default = 10.0, metavar = 'PCT',
        help = 'Percentile level for clipping color scale.')

    args = parser.parse_args()
    if args.no_display and not args.output_name:
        print 'No display our output requested.'
        return 0
    if args.galaxy is None and args.group is None:
        print 'Must specify either a galaxy or a group.'
        return -1
    if args.galaxy is not None and args.group is not None:
        print 'Cannot specify both a galaxy and a group.'
        return -1
    if args.partials + args.matrix + args.covariance + args.correlation > 1:
        print 'Can only specify one of the partials,matrix,covariance options.'
        return -1

    if args.clip_percentile < 0 or args.clip_percentile >= 50:
        print 'Invalid --clip-percentile %f (should be 0-50).' % args.clip_percentile
        return -1

    # Load the analysis results file we will get partial derivative images from.
    try:
        reader = descwl.output.Reader.from_args(defer_stamp_loading = True,args = args)
        results = reader.results
        npartials = len(results.slice_labels)
        if args.verbose:
            print results.survey.description()
    except RuntimeError,e:
        print str(e)
        return -1
    if results.table is None:
        print 'Input file is missing a results catalog.'
        return -1
    if results.stamps is None:
        print 'Input file is missing stamp datacubes.'
        return -1

    # Look for the selected galaxy or group.
    if args.galaxy:
        selected = results.select('db_id==%d' % args.galaxy)
        if len(selected) == 0:
            print 'No such galaxy with ID %d' % args.galaxy
            return -1
        title = 'galaxy-%d' % args.galaxy
    else:
        selected = results.select('grp_id==%d' % args.group)
        if len(selected) == 0:
            print 'No such group with ID %d' % args.group
            return -1
        title = 'group-%d' % args.group
    # Sort selected galaxies in increasing rank order.
    sort_order = np.argsort(results.table['grp_rank'][selected])
    selected = selected[sort_order]
    num_selected = len(selected)
    npar = npartials*num_selected
    nrows,ncols = npar,npar

    # Get the background image for these galaxies.
    background = results.get_subimage(selected)
    height,width = background.array.shape

    # Calculate matrix elements.
    fisher,covariance,variance,correlation = results.get_matrices(selected)
    show_matrix = args.matrix or args.covariance or args.correlation
    if show_matrix:
        if args.matrix:
            matrix = fisher
        elif args.covariance:
            matrix = covariance
        else:
            matrix = correlation

    #do bias calculations
    bias = results.get_bias(selected, covariance)

    if args.bias:
        import math
        #print out galaxy id + each of the params
        slice_labels = ['flux','x','y','s','g1','g2']
        for i in range(len(bias)):
            print slice_labels[i]
            print 'std:', math.sqrt(covariance[i][i])
            print 'bias:', bias[i]
            print 'bias/std:',bias[i]/math.sqrt(covariance[i][i])

        # print selected 
        # print covariance
        # print bias

    # Print a summary table of RMS errors on each parameter.
    if args.verbose and correlation is not None:
        dtypes = [ (name,np.float32) for name in results.slice_labels ]
        dtypes.insert(0,('id',np.int64))
        summary = np.empty(shape = (len(selected),),dtype = dtypes)
        summary['id'] = results.table['db_id'][selected]
        for index in range(ncols):
            galaxy = index//npartials
            islice = index%npartials
            summary[galaxy][islice+1] = np.sqrt(variance[index])
        print astropy.table.Table(summary)

    # Calculate the bounds for our figure.
    if args.partials:
        nrows = 1
    figure_scale = args.figure_size/(ncols*max(height,width))
    figure_width = ncols*width*figure_scale
    figure_height = nrows*height*figure_scale
    figure = plt.figure(figsize = (figure_width,figure_height),frameon=False)
    figure.canvas.set_window_title(title)
    plt.subplots_adjust(left = 0,bottom = 0,right = 1,top = 1,wspace = 0,hspace = 0)

    def draw(row,col,pixels):
        axes = plt.subplot(nrows,ncols,row*ncols+col+1)
        axes.set_axis_off()
        if row == col:
            # All values are positive.
            vmax = np.percentile(pixels[pixels != 0], 100 - args.clip_percentile)
        else:
            vmax = np.max(np.fabs(np.percentile(pixels[pixels != 0],
                (args.clip_percentile, 100 - args.clip_percentile))))
        vmin = -vmax
        scaled = np.clip(pixels,vmin,vmax)
        plt.imshow(scaled,origin = 'lower',interpolation = 'nearest',
            cmap = args.colormap,vmin = vmin,vmax = vmax)

    def draw_param_label(index,row,col):
        # index determines which parameter label to draw.
        # row,col determine where the label will be drawn.
        islice = index%npartials
        igalaxy = index//npartials
        rank = results.table['grp_rank'][selected[igalaxy]]
        # Latex labels do not get the correct vertical alignment.
        ##param_label = '$%s_{%d}$' % (results.slice_labels[islice],rank)
        param_label = '%s_%d' % (results.slice_labels[islice],rank)
        x = (col+1.)/ncols
        y = 1.-float(row)/nrows
        plt.annotate(param_label,xy = (x,y),xycoords = 'figure fraction',
        color = args.label_color, fontsize = args.label_size,
        horizontalalignment = 'right',verticalalignment = 'top')

    if args.partials:
        # Draw the partial-derivative images on a single row.
        stamp = results.get_subimage(selected)
        for col in range(ncols):
            galaxy = selected[col//npartials]
            islice = col%npartials
            stamp.array[:] = 0.
            stamp[results.bounds[galaxy]] = results.get_stamp(galaxy,islice)
            if islice == 0:
                # Normalize to give partial with respect to added flux in electrons.
                stamp /= results.table['flux'][galaxy]
            draw(0,col,stamp.array)
            if not args.no_labels:
                draw_param_label(index=col,row=0,col=col)
    elif show_matrix:
        # Draw the values of the matrix we calculated above.
        span = np.arange(nrows)
        row,col = np.meshgrid(span,span)
        lower_triangle = np.ma.masked_where(row > col,matrix)
        axes = plt.subplot(1,1,1)
        axes.set_axis_off()
        vmin,vmax = (-1.,+1.) if args.correlation else (None,None)
        plt.imshow(lower_triangle,interpolation = 'nearest',aspect = 'auto',
            cmap = args.colormap,vmin = vmin,vmax = vmax)
        if not args.no_labels:
            for row in range(nrows):
                for col in range(row+1):
                    value_label = args.value_format % matrix[row,col]
                    xc = (col+0.5)/ncols
                    yc = 1.-(row+0.5)/nrows
                    plt.annotate(value_label,xy = (xc,yc),xycoords = 'figure fraction',
                        color = args.label_color, fontsize = args.label_size,
                        horizontalalignment = 'center',verticalalignment = 'center')
                    if row == col and not args.no_labels:
                        draw_param_label(index=row,row=row,col=col)
    else:
        # Draw Fisher-matrix images.
        stamp = background.copy()
        for row,index1 in enumerate(selected):
            for col,index2 in enumerate(selected[:row+1]):
                images,overlap = results.get_fisher_images(index1,index2,background)
                if overlap is None:
                    continue
                for par1 in range(npartials):
                    fisher_row = npartials*row+par1
                    for par2 in range(npartials):
                        fisher_col = npartials*col+par2
                        if fisher_col > fisher_row:
                            continue
                        stamp.array[:] = 0.
                        stamp[overlap].array[:] = images[par1,par2]
                        draw(fisher_row,fisher_col,stamp.array)
                        if fisher_row == fisher_col and not args.no_labels:
                            draw_param_label(index = fisher_row,row = fisher_row,col = fisher_col)

    if args.output_name:
        figure.savefig(args.output_name)

    if not args.no_display:
        plt.show()

if __name__ == '__main__':
    main()
