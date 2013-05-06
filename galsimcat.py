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

# image pixel bounds
xmin = ymin = 1792
xmax = ymax = 2304
##xmin = ymin = 0
##xmax = ymax = 50

pixelScale = 0.2 # arcsecs / pixel

# Calculate RA,DEC bounds
RAmin = xmin*pixelScale
RAmax = xmax*pixelScale
DECmin = ymin*pixelScale
DECmax = ymax*pixelScale

# Number of visits to simulate (source flux will be multiplied by 2*nvisits)
nvisits = 230

# average sky level to simulate in ADU units per pixel
# see http://darkmatter.ps.uci.edu/wiki/LSST%20Exposure%20Calculations
skyLevel = 780.778

# Specify constant cosmic shear to apply across the field
g1 = 0.0
g2 = 0.0

# Configure parameter variations
dxyc = pixelScale/3.
nxyc = 1
dg12 = 0.03
ng12 = 1
dhlr = 0.02
nhlr = 1

def createSource(flux,rhalf,q,beta,g1,g2,dx,dy):
    source = galsim.Exponential(flux = flux, half_light_radius = rhalf)
    ##source = galsim.Gaussian(flux = flux, half_light_radius = rhalf)
    source.applyShear(q = q, beta = beta*galsim.radians)
    source.applyShear(g1 = g1, g2 = g2)
    source.applyShift(dx = dx, dy = dy)
    return source

def main(argv):
    
    # Parse command-line args
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", default='trim.dat', help="name of input catalog to read")
    args = parser.parse_args()

    # In non-script code, use getLogger(__name__) at module scope instead.
    logging.basicConfig(format="%(message)s", level=logging.INFO, stream=sys.stdout)
    logger = logging.getLogger("aliasing")

    # Define the image characteristics
    pix = galsim.Pixel(pixelScale)

    # Open the input catalog to use
    cat = galsim.InputCatalog(args.input)
    logger.info('Reading input catalog %r' % args.input)
    
    # Initialize the list of stamps we will create and an empty image that
    # we will add each stamp to
    stamps = [ ]
    field = galsim.ImageD(xmax-xmin,ymax-ymin,init_value = 0)

    # Define the psf (beta = 3 taken from GREAT10 simulations)
    psf = galsim.Moffat(beta = 3, fwhm = 0.7)

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
        flux = 2*nvisits*flux
        
        # Calculate the offset of this stamp's central pixel (in pixels) in the full field
        xc = int(round(RA/pixelScale))
        yc = int(round(DEC/pixelScale))

        # Calculate the subpixel shift of the object's center (in arcsecs) within the stamp
        xshift = RA-pixelScale*xc
        yshift = DEC-pixelScale*yc

        # Calculate the stamp size to use. Always round up to an odd integer
        # so that flux is consistently centered (see Issue #380).
        w = 2*int(math.ceil(width/pixelScale))+1
        h = 2*int(math.ceil(height/pixelScale))+1
        ## make the stamp square
        #w = h = max(w,h)
        logger.info('rendering stamp %d (id %d) with size %d x %d' % (nkeep,k,w,h))

        # Calculate the amount to shift this stamp so that it is correctly centered
        # in our image with bounds (xmin,ymin) - (xmax,ymax).
        ## Do these need an addition (1,1) offset?
        dx = xc - (w-1)/2 - xmin
        dy = yc - (h-1)/2 - ymin

        # Define the instrinsic galaxy profile
        gal = createSource(flux,hlr,q,beta,g1,g2,xshift,yshift)
        
        # Create a copy that will not have any psf applied
        nopsf = createSource(flux,hlr,q,beta,g1,g2,xshift,yshift)

        # Create a family of galaxies with small parameter variations
        variations = [ nopsf, gal ]
        for ixyc in range(1,nxyc+1):
            delta = ixyc*dxyc
            # Vary xc by +delta
            copy = createSource(flux,hlr,q,beta,g1,g2,xshift+delta,yshift)
            variations.append(copy)
            # Vary xc by -delta
            copy = createSource(flux,hlr,q,beta,g1,g2,xshift-delta,yshift)
            variations.append(copy)
            # Vary yc by +delta
            copy = createSource(flux,hlr,q,beta,g1,g2,xshift,yshift+delta)
            variations.append(copy)
            # Vary yc by -delta
            copy = createSource(flux,hlr,q,beta,g1,g2,xshift,yshift-delta)
            variations.append(copy)
        for ie12 in range(1,ng12+1):
            delta = ie12*dg12
            # Vary e1 by +delta
            copy = createSource(flux,hlr,q,beta,g1+delta,g2,xshift,yshift)
            variations.append(copy)
            # Vary e1 by -delta
            copy = createSource(flux,hlr,q,beta,g1-delta,g2,xshift,yshift)
            variations.append(copy)
            # Vary e2 by +delta
            copy = createSource(flux,hlr,q,beta,g1,g2+delta,xshift,yshift)
            variations.append(copy)
            # Vary e2 by -delta
            copy = createSource(flux,hlr,q,beta,g1,g2-delta,xshift,yshift)
            variations.append(copy)
        for ihlr in range(1,nhlr+1):
            delta = ihlr*dhlr
            # Vary hlr by +delta
            copy = createSource(flux,hlr+delta,q,beta,g1,g2,xshift,yshift)
            variations.append(copy)
            # Vary hlr by -delta
            copy = createSource(flux,hlr-delta,q,beta,g1,g2,xshift,yshift)
            variations.append(copy)

        # Loop over variations to render
        for src in variations:
            stamp = galsim.ImageD(w,h)
            if src == nopsf:
                object = galsim.Convolve([src,pix], real_space = True)
            else:
                object = galsim.Convolve([src,psf,pix])
            object.draw(image = stamp, dx = pixelScale)
            stamp.shift(dx,dy)
            stamps.append(stamp)
            if src == gal:
                # Add the nominal galaxy to the full field image
                overlap = stamp.bounds & field.bounds
                field[overlap] += stamp[overlap]

    # Write the full field image to a separate file
    outname = 'field.fits'
    logger.info('Saving full field to %r' % outname)
    galsim.fits.write(field,outname)

    # Write out the full field image with noise added
    if skyLevel > 0:
        rng = galsim.BaseDeviate(123)
        noise = galsim.PoissonNoise(rng,sky_level = 2*nvisits*skyLevel)
        field.addNoise(noise)
        outname = 'fieldnoise.fits'
        logger.info('Saving full field to %r' % outname)
        galsim.fits.write(field,outname)

    # Write the indiviual stamps to multifits files
    outname = 'each.fits'
    logger.info('Saving %d source stamps to %r' % (len(stamps),outname))
    galsim.fits.writeMulti(stamps, outname)

if __name__ == "__main__":
    main(sys.argv)
