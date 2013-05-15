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

"""
Creates a source object with the specified parameters.
"""
def createSource(flux,xc,yc,hlr,q,beta,g1,g2):
    source = galsim.Exponential(flux = flux, half_light_radius = hlr)
    source.applyShear(q = q, beta = beta*galsim.radians)
    source.applyShear(g1 = g1, g2 = g2)
    source.applyShift(dx = xc, dy = yc)
    return source

"""
Renders the specified source convolved with a psf (which might be None)
and pixel response into a postage stamp with the specified bounding box.
"""
def createStamp(src,psf,pix,bbox):
    stamp = galsim.ImageD(bbox)
    if psf == None:
        obj = galsim.Convolve([src,pix], real_space = True)
    else:
        obj = galsim.Convolve([src,psf,pix])
    obj.draw(image = stamp)
    return stamp

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
    parser.add_argument("--psf-fwhm", type = float, default = 0.7,
        help = "psf full-width-half-max in arcsecs")
    parser.add_argument("--psf-beta", type = float, default = 3.0,
        help = "psf Moffat parameter beta")
    parser.add_argument("--nvisits", type = int, default = 230,
        help = "number of visits to simulate")
    parser.add_argument("--sky-level", type = float, default = 780.778,
        help = "average sky level to simulate (ADU/pixel)")
    parser.add_argument("--g1", type = float, default = 0.,
        help = "Constant shear component g1 to apply")
    parser.add_argument("--g2", type = float, default = 0.,
        help = "Constant shear component g2 to apply")
    parser.add_argument("--partials", action = "store_true",
        help = "Calculate and save partial derivatives with respect to object parameters")
    parser.add_argument("--partials-order", type = int, default = 1,
        help = "Order of finite difference equation to use for evaluating partials")
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

    # Initialize finite difference calculations if necessary
    if args.partials:
        if args.partials_order < 1 or args.partials_order > 4:
            logger.error('Bad parameter: partials-order must be an integer 1-4.')
            sys.exit(-1)
        # Initialize the finite difference coefficients to use
        if args.partials_order == 1:
            fdCoefs = (1./2.,)
        elif args.partials_order == 2:
            fdCoefs = (2./3.,-1./12.)
        elif args.partials_order == 3:
            fdCoefs = (3./4.,-3./20.,1./60.)
        else:
            fdCoefs = (4./5.,-1./5.,4./105.,-1./280.)

    # Define the image characteristics
    pix = galsim.Pixel(args.pixel_scale)

    # Open the input catalog to use
    cat = galsim.InputCatalog(args.input)
    logger.info('Reading input catalog %r' % args.input)
    
    # Define the psf to use
    psf = galsim.Moffat(beta = args.psf_beta, fwhm = args.psf_fwhm)

    # Create an empty image that we will add each stamp to
    field = galsim.ImageD(args.xmax-args.xmin,args.ymax-args.ymin,init_value = 0)

    # Initialize the list of per-object stamp HDUs we will fill
    hdu = pyfits.PrimaryHDU()
    hduList = pyfits.HDUList([hdu])

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
        logger.info('rendering stamp %d (id %d) with size %d x %d' % (nkeep,k,w,h))
        
        # Calculate the amount to shift this stamp so that it is correctly centered
        # in our image with bounds (xmin,ymin) - (xmax,ymax).
        ## Do these need an addition (1,1) offset?
        dx = xc - (w-1)/2 - args.xmin
        dy = yc - (h-1)/2 - args.ymin
        
        # Combined the stamp size and offset into a bounding box for this stamp
        bbox = galsim.BoundsI(dx,dx+w-1,dy,dy+h-1)

        # Define the nominal source parameters for this object
        params = {
            'flux':flux, 'xc':xshift, 'yc':yshift,
            'hlr':hlr, 'q':q, 'beta':beta,
            'g1':args.g1, 'g2': args.g2
        }

        # Create stamps for the galaxy with and w/o the psf applied
        gal = createSource(**params)
        nopsf = createStamp(gal,None,pix,bbox)
        nominal = createStamp(gal,psf,pix,bbox)

        # Add the nominal galaxy to the full field image
        overlap = nominal.bounds & field.bounds
        field[overlap] += nominal[overlap]
        
        # Initialize the datacube of stamps that we will save for this object
        datacube = [ nopsf, nominal ]

        if args.partials:
            # Specify the amount to vary each parameter for partial derivatives
            # (we don't use a dictionary here since we want to control the order)
            variations = [
                ('xc',args.pixel_scale/3.), ('yc',args.pixel_scale/3.),
                ('hlr',0.05*hlr),
                ('g1',0.03), ('g2',0.03)
            ]
            # loop over source parameters to vary
            for (pname,delta) in variations:
                # create stamps for each variation of this parameter
                newparams = params.copy()
                partial = galsim.ImageD(w,h)
                for step in range(args.partials_order):
                    # create and save the positive variation stamp
                    newparams[pname] = params[pname] + (step+1)*delta
                    newsource = createSource(**newparams)
                    plus = createStamp(newsource,psf,pix,bbox)
                    # create and save the negative variation stamp
                    newparams[pname] = params[pname] - (step+1)*delta
                    newsource = createSource(**newparams)
                    minus = createStamp(newsource,psf,pix,bbox)
                    # update the finite difference calculation of this partial
                    partial += (fdCoefs[step]/delta)*(plus - minus)
                # append this partial to our datacube
                datacube.append(partial)

        # Add a new HDU with a datacube for this object's stamps
        # We don't use compression = 'gzip_tile' for now since it is lossy
        # and mathematica cannot Import it.
        galsim.fits.writeCube(datacube, hdu_list = hduList)

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
    galsim.fits.write_file(outname, hdus = hduList, clobber = True,
        file_compress = None, pyfits_compress = None)

if __name__ == "__main__":
    main()
