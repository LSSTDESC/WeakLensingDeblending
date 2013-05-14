#!/usr/bin/env python
#######################################################################################
## TODO:
## - add bulge components
## - use more realistic (COSMOS) distribution of intrinsic ellipticities?
## - add variations over galaxy size parameters (and others?)
#######################################################################################

import sys
import os
import math
import argparse

import logging
import galsim
import pyfits

def createSource(flux,rhalf,q,beta,g1,g2,dx,dy):
    source = galsim.Exponential(flux = flux, half_light_radius = rhalf)
    ##source = galsim.Gaussian(flux = flux, half_light_radius = rhalf)
    source.applyShear(q = q, beta = beta*galsim.radians)
    source.applyShear(g1 = g1, g2 = g2)
    source.applyShift(dx = dx, dy = dy)
    return source

def main():
   
    # Parse command-line args
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", default = 'trim.dat',
        help = "name of input catalog to read")
    parser.add_argument("-o","--output", default = 'catout',
        help = "base name of output files to write")
    parser.add_argument("--xmin", type = int, default = 1792,
        help = "left edge of image (pixels)")
    parser.add_argument("--xmax", type = int, default = 2304,
        help = "right edge of image (pixels)")
    parser.add_argument("--ymin", type = int, default = 1792,
        help = "bottom edge of image (pixels)")
    parser.add_argument("--ymax", type = int, default = 2304,
        help = "top edge of image (pixels)")
    parser.add_argument("--pixel-scale", type = float, default = 0.2,
        help = "pixel scale (arscecs/pixel)")
    parser.add_argument("--nvisits", type = int, default = 230,
        help = "number of visits to simulate")
    parser.add_argument("--sky-level", type = float, default = 780.778,
        help = "average sky level to simulate (ADU/pixel)")
    parser.add_argument("--g1", type = float, default = 0.,
        help = "Constant shear component g1 to apply")
    parser.add_argument("--g2", type = float, default = 0.,
        help = "Constant shear component g2 to apply")
    args = parser.parse_args()

    # In non-script code, use getLogger(__name__) at module scope instead.
    logging.basicConfig(format="%(message)s", level=logging.INFO, stream=sys.stdout)
    logger = logging.getLogger("aliasing")
    logger.info('Using output prefix %r' % args.output)

    # Calculate RA,DEC bounds
    RAmin = args.xmin*args.pixel_scale
    RAmax = args.xmax*args.pixel_scale
    DECmin = args.ymin*args.pixel_scale
    DECmax = args.ymax*args.pixel_scale

    # Configure parameter variations
    dxyc = args.pixel_scale/3.
    nxyc = 1
    dg12 = 0.03
    ng12 = 1
    dhlr = 0.02
    nhlr = 1

    # Define the image characteristics
    pix = galsim.Pixel(args.pixel_scale)

    # Open the input catalog to use
    cat = galsim.InputCatalog(args.input)
    logger.info('Reading input catalog %r' % args.input)
    
    # Define the psf (beta = 3 taken from GREAT10 simulations)
    psf = galsim.Moffat(beta = 3, fwhm = 0.7)

    # Create an empty image that we will add each stamp to
    field = galsim.ImageD(args.xmax-args.xmin,args.ymax-args.ymin,init_value = 0)

    # Initialize the list of per-object stamp HDUs we will fill
    hdu = pyfits.PrimaryHDU()
    stampsList = pyfits.HDUList([hdu])

    # Loop over catalog entries
    nkeep = 0
    for k in range(cat.nobjects):

        RA = cat.getFloat(k,0)     # arcsecs
        DEC = cat.getFloat(k,1)    # arcsecs
        flux = cat.getFloat(k,2)   # integrated ADU / exposure
        hlr = cat.getFloat(k,3)    # arcsecs
        q = cat.getFloat(k,4)      # 0 < q < 1
        beta = cat.getFloat(k,5)   # degrees
        width = cat.getFloat(k,6)  # half-width of bounding box in arcsecs
        height = cat.getFloat(k,7) # half-height of bounding box in arcsecs
        iAB = cat.getFloat(k,8)    # AB mags
        
        # Is this object's center within our image?
        if RA < RAmin or RA > RAmax or DEC < DECmin or DEC > DECmax:
            continue
        nkeep += 1
        
        # Scale flux to number of vists (extra factor of 2 because 1 visit = 2 exposures)
        flux = 2*args.nvisits*flux
        
        # Calculate the offset of this stamp's central pixel (in pixels) in the full field
        xc = int(round(RA/args.pixel_scale))
        yc = int(round(DEC/args.pixel_scale))

        # Calculate the subpixel shift of the object's center (in arcsecs) within the stamp
        xshift = RA-args.pixel_scale*xc
        yshift = DEC-args.pixel_scale*yc

        # Calculate the stamp size to use. Always round up to an odd integer
        # so that flux is consistently centered (see Issue #380).
        w = 2*int(math.ceil(width/args.pixel_scale))+1
        h = 2*int(math.ceil(height/args.pixel_scale))+1
        ## make the stamp square
        #w = h = max(w,h)
        logger.info('rendering stamp %d (id %d) with size %d x %d' % (nkeep,k,w,h))

        # Calculate the amount to shift this stamp so that it is correctly centered
        # in our image with bounds (xmin,ymin) - (xmax,ymax).
        ## Do these need an addition (1,1) offset?
        dx = xc - (w-1)/2 - args.xmin
        dy = yc - (h-1)/2 - args.ymin

        # Define the instrinsic galaxy profile
        gal = createSource(flux,hlr,q,beta,args.g1,args.g2,xshift,yshift)
        
        # Create a copy that will not have any psf applied
        nopsf = createSource(flux,hlr,q,beta,args.g1,args.g2,xshift,yshift)

        # Create a family of galaxies with small parameter variations
        variations = [ nopsf, gal ]
        for ixyc in range(1,nxyc+1):
            delta = ixyc*dxyc
            # Vary xc by +delta
            copy = createSource(flux,hlr,q,beta,args.g1,args.g2,xshift+delta,yshift)
            variations.append(copy)
            # Vary xc by -delta
            copy = createSource(flux,hlr,q,beta,args.g1,args.g2,xshift-delta,yshift)
            variations.append(copy)
            # Vary yc by +delta
            copy = createSource(flux,hlr,q,beta,args.g1,args.g2,xshift,yshift+delta)
            variations.append(copy)
            # Vary yc by -delta
            copy = createSource(flux,hlr,q,beta,args.g1,args.g2,xshift,yshift-delta)
            variations.append(copy)
        for ie12 in range(1,ng12+1):
            delta = ie12*dg12
            # Vary e1 by +delta
            copy = createSource(flux,hlr,q,beta,args.g1+delta,args.g2,xshift,yshift)
            variations.append(copy)
            # Vary e1 by -delta
            copy = createSource(flux,hlr,q,beta,args.g1-delta,args.g2,xshift,yshift)
            variations.append(copy)
            # Vary e2 by +delta
            copy = createSource(flux,hlr,q,beta,args.g1,args.g2+delta,xshift,yshift)
            variations.append(copy)
            # Vary e2 by -delta
            copy = createSource(flux,hlr,q,beta,args.g1,args.g2-delta,xshift,yshift)
            variations.append(copy)
        for ihlr in range(1,nhlr+1):
            delta = ihlr*dhlr
            # Vary hlr by +delta
            copy = createSource(flux,hlr+delta,q,beta,args.g1,args.g2,xshift,yshift)
            variations.append(copy)
            # Vary hlr by -delta
            copy = createSource(flux,hlr-delta,q,beta,args.g1,args.g2,xshift,yshift)
            variations.append(copy)

        # Loop over variations to render
        cubeStamps = [ ]
        for src in variations:
            stamp = galsim.ImageD(w,h)
            if src == nopsf:
                object = galsim.Convolve([src,pix], real_space = True)
            else:
                object = galsim.Convolve([src,psf,pix])
            object.draw(image = stamp, dx = args.pixel_scale)
            stamp.shift(dx,dy)
            cubeStamps.append(stamp)
            if src == gal:
                # Add the nominal galaxy to the full field image
                overlap = stamp.bounds & field.bounds
                field[overlap] += stamp[overlap]
                
        # Add a new HDU with a datacube for this object's stamps
        # We don't use compression = 'gzip_tile' for now since it is lossy
        # and mathematica cannot Import it.
        galsim.fits.writeCube(cubeStamps, hdu_list = stampsList)

    # Write the full field image to a separate file
    outname = args.output + '_field.fits'
    logger.info('Saving full field to %r' % outname)
    galsim.fits.write(field,outname)

    # Write out the full field image with noise added
    if args.sky_level > 0:
        rng = galsim.BaseDeviate(123)
        noise = galsim.PoissonNoise(rng,sky_level = 2*args.nvisits*args.sky_level)
        field.addNoise(noise)
        outname = args.output + '_noise.fits'
        logger.info('Saving full field to %r' % outname)
        galsim.fits.write(field,outname)

    # Write the object stamp datacubes
    outname = args.output + '_stamps.fits'
    logger.info('Saving stamps to %r' % outname)
    galsim.fits.write_file(outname, hdus = stampsList, clobber = True,
        file_compress = None, pyfits_compress = None)

if __name__ == "__main__":
    main()
