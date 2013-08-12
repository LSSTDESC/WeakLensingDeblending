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

twopi = 2*math.pi
deg2rad = math.pi/180.
deg2arcsec = 3600.
deg2arcmin = 60.

"""
Creates a source object with the specified parameters.
"""
def createSource(flux,bulgeFraction,xc,yc,hlr_d,q_d,beta_d,hlr_b,q_b,beta_b,g1,g2):
    # Define the disk component, if any
    if bulgeFraction < 1:
        disk = galsim.Exponential(flux = flux*(1-bulgeFraction), half_light_radius = hlr_d)
        disk.applyShear(q = q_d, beta = beta_d*galsim.radians)
    # Define the bulge component, if any
    if bulgeFraction > 0:
        bulge = galsim.DeVaucouleurs(flux = flux*bulgeFraction, half_light_radius = hlr_b)
        bulge.applyShear(q = q_b, beta = beta_b*galsim.radians)
    # Combined the disk and bulge components
    if bulgeFraction == 0:
        source = disk
    elif bulgeFraction == 1:
        source = bulge
    else:
        source = disk + bulge
    # Shift and shear the combined source
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
        gsp = galsim.GSParams(maximum_fft_size=8192)
        obj = galsim.Convolve([src,psf,pix],gsparams=gsp)
    obj.draw(image = stamp, dx = pix.getXWidth())
    return stamp

"""
Returns (dx,dy) for the bounding box of a surface brightness profile
SB(x,y) whose isophotes are ellipses with the shape (q,beta) and which
has an underlying normalized radial profile p(r). The inputs are:

  maxSB = totalFlux*p(0) = maximum surface brightness before shear
  thresholdSB = threshold surface brightness after shear
  q = ratio of minor to major axes of ellipse with 0 < q <= 1
  beta = angle of ellipse's major axis in radians
  rFunction = returns R(b) such that p(R) = b*p(0) with 0 < b < 1

The returned (dx,dy) are in arcsecs, and defined such that SB(x,y) < f0
is guaranteed for |x| > dx or |y| > dy. The result only depends on the
ratio thresholdSB/maxSB so they must be in the same (abitrary) units.
"""
def boundingBox(maxSB,thresholdSB,q,beta,rFunction):
    # Calculate shear affine transform parameters
    g = (1-q)/(1+q)
    gp = g*math.cos(2*beta)
    gx = g*math.sin(2*beta)
    detM = 1 - gp*gp - gx*gx
    # Calculate the dimensionless surface brightness ratio at threshold.
    b = thresholdSB/(maxSB*detM)
    if b <= 0:
        raise RuntimeError('boundingBox: invalid input parameters')
    if b >= 1:
        # The max surface brightness is below our threshold SB(0,0) <= f0
        return (0,0)
    # Calculate the threshold radius of the radial profile.
    rcut = rFunction(b)
    # Shear this circle and return its bounding box dimensions
    dx = rcut*math.sqrt(((1+gp)*(1+gp)+gx*gx)/detM) # half width in arcsecs
    dy = rcut*math.sqrt(((1-gp)*(1-gp)+gx*gx)/detM) # half height in arcsecs
    return (dx,dy)

"""
Returns (dx,dy) for the bounding box of a Sersic profile with n = 1 or 4.
The input flux should be in ADU, hlr in arscecs, beta in radians, f0 in
ADU/arcsec^2. 0 < q <= 1 is dimensionless. The returned (dx,dy) are in
arcsecs. See boundingBox above for details.
"""
def sersicBounds(n,flux,hlr,q,beta,f0):
    # Convert the half-light radius to the appropriate scale radius r0
    # and calculate the Sersic normalization constant
    if n == 1:
        r0 = hlr/1.67835
        norm = twopi*r0*r0
    elif n == 4:
        r0 = hlr/3459.49
        norm = 20160*twopi*r0*r0  # 20160 = n*Gamma[2*n]
    else:
        raise RuntimeError('Sersic index n = %d is not supported.' % n)
    # Calculate and return the bounding box
    return boundingBox(flux/norm,f0,q,beta,
        lambda b: r0*math.pow(-math.log(b),n))

"""
Returns (dx,dy) for the bounding box of a Moffat profile. The input flux
should be in ADU, fwhm in arcsecs, beta in radians, f0 in ADU/arcsec^2.
0 < q <= 1 and beta > 1 are dimensionless. The returned (dx,dy) are in
arcsecs. See boundingBox above for details.
"""
def moffatBounds(moffatBeta,flux,fwhm,q,beta,f0):
    # Check that beta is valid
    if moffatBeta <= 1:
        raise RuntimeError('Moffat beta < 1 is not valid.')
    # Convert the fwhm to the corresponding scale radius
    r0 = 0.5*fwhm/math.sqrt(math.pow(2,1./moffatBeta)-1)
    # Calculate the normalization factor norm = 1/p(0)
    norm = math.pi*r0*r0/(moffatBeta-1)
    # Calculate and return the bounding box
    return boundingBox(flux/norm,f0,q,beta,
        lambda b: r0*math.sqrt(1-math.pow(b,(moffatBeta-1.)/moffatBeta)))
    
"""
Returns a mask image of values 0 or 1 depending on whether the corresponding
input image pixel value is above or below the specified threshold in ADU.
"""
def createMask(image,threshold,args):
    # create an empty mask image with the same dimensions as the input image
    box = image.bounds
    mask = galsim.ImageS(box)
    if not args.no_trim:
        # initialize our trimmed bounds to just the central pixel
        # (the numerator should always be even for odd width,height)
        xmin = (box.getXMin()+box.getXMax())
        ymin = (box.getYMin()+box.getYMax())
        assert xmin%2 == 0 and ymin%2 == 0
        xmin = xmin/2
        ymin = ymin/2
        xmax = xmin
        ymax = ymin
    # loop over image pixels
    for (rowIndex,row) in enumerate(image.array):
        y = box.getYMin()+rowIndex
        for (pixelIndex,pixelValue) in enumerate(row):
            x = box.getXMin()+pixelIndex
            if pixelValue >= threshold:
                mask.array[rowIndex,pixelIndex] = 1
                if not args.no_trim:
                    xmin = min(x,xmin)
                    xmax = max(x,xmax)
                    ymin = min(y,ymin)
                    ymax = max(y,ymax)
    if not args.no_trim:
        trimmed = galsim.BoundsI(xmin,xmax,ymin,ymax)
        mask = mask[trimmed]
    return mask

"""
Performs any final processing on stamp, controlled by args, then appends it to stamps.
Saved stamps are always normalized to a single exposure, i.e., they are scaled by
1/(2*nvisits) relative to what gets added to the full field image.
Returns True if the stamp was saved, or otherwise False.
"""
def saveStamp(stamps,stamp,trimmed,args):
    # Trim the stamp to its threshold bounding box
    if not args.no_trim:
        stamp = stamp[trimmed]
    # Clip the stamp so that does not extend beyond the field image. This results
    # in potentially smaller files with sources that might not be centered.
    if not args.no_clip:
        overlap = stamp.bounds & galsim.BoundsI(1,args.width,1,args.height)
        if overlap.area() == 0:
            # skip this stamp if it falls completely outside our field (after trimming)
            return False
        stamp = stamp[overlap]
    # Remember this stamp.
    stamps.append(stamp)
    return True

def main():

    # Parse command-line args
    parser = argparse.ArgumentParser()
    parser.add_argument("-v", "--verbose", action = "store_true",
        help = "provide more verbose output")
    parser.add_argument("-i", "--input", default = 'gcat.dat',
        help = "name of input catalog to read")
    parser.add_argument("-o","--output", default = 'catout',
        help = "base name of output files to write")
    parser.add_argument("-x","--x-center", type = float, default = 0.5,
        help = "central RA of image (degrees)")
    parser.add_argument("-y","--y-center", type = float, default = 0.0,
        help = "central DEC of image (degrees)")
    parser.add_argument("--width", type = int, default = 512,
        help = "image width (pixels)")
    parser.add_argument("--height", type = int, default = 512,
        help = "image height (pixels)")
    parser.add_argument("--margin", type = float, default = 10.,
        help = "size of surrounding margin where objects might leak into image (arcmins)")
    parser.add_argument("--pixel-scale", type = float, default = 0.2,
        help = "pixel scale (arscecs/pixel)")
    parser.add_argument("--psf-fwhm", type = float, default = 0.7,
        help = "psf full-width-half-max in arcsecs (zero for no psf)")
    parser.add_argument("--psf-beta", type = float, default = 3.0,
        help = "psf Moffat parameter beta")
    parser.add_argument("--band", choices = ['u','g','r','i','z','y'], default = 'i',
        help = "LSST imaging band to use for source fluxes")
    parser.add_argument("--flux-norm", type = float, default = 711.,
        help = "total flux in ADU for one exposure of a typical galaxy of AB mag 24")
    parser.add_argument("--sky-level", type = float, default = 780.778,
        help = "average sky level to simulate (ADU/pixel in one exposure)")
    parser.add_argument("--sn-cut", type = float, default = 0.1,
        help = "keep all pixels above this signal-to-noise ratio cut")
    parser.add_argument("--nvisits", type = int, default = 230,
        help = "number of visits to simulate (1 visit = 2 exposures)")
    parser.add_argument("--g1", type = float, default = 0.,
        help = "constant shear component g1 to apply")
    parser.add_argument("--g2", type = float, default = 0.,
        help = "constant shear component g2 to apply")
    parser.add_argument("--stamps", action = "store_true",
        help = "save postage stamps for each source (normalized to 1 exposure)")
    parser.add_argument("--no-clip", action = "store_true",
        help = "do not clip stamps to the image bounds")
    parser.add_argument("--no-trim", action = "store_true",
        help = "do not trim stamps to their threshold bounding box")
    parser.add_argument("--no-bulge", action = "store_true",
        help = "do not include any galactic bulge components")
    parser.add_argument("--partials", action = "store_true",
        help = "calculate and save partial derivatives wrt shape parameters (normalized to 1 exposure)")
    parser.add_argument("--partials-order", type = int, default = 1,
        help = "order of finite difference equation to use for evaluating partials")
    parser.add_argument("--only-line", type = int, default = 0,
        help = "only use the specified line number from the input catalog (when non-zero)")
    args = parser.parse_args()

    # Configure the GalSim logger
    logging.basicConfig(format="%(message)s", level=logging.INFO, stream=sys.stdout)
    logger = logging.getLogger("galsimcat")
    logger.info('Using output prefix %r' % args.output)

    # Define the pixel response
    pix = galsim.Pixel(args.pixel_scale)

    # Define the psf to use
    if args.psf_fwhm > 0:
        psf = galsim.Moffat(beta = args.psf_beta, fwhm = args.psf_fwhm)
    else:
        psf = None

    # Create an empty image that represents the whole field
    field = galsim.ImageD(args.width,args.height)
    field.setScale(pix.getXWidth())
    
    # Calculate the relative scaling of RA and angles relative to the image center
    RAscale = math.cos(args.y_center*deg2rad)

    # Calculate the corners of the image in degrees
    RAmin = args.x_center - 0.5*args.width*args.pixel_scale/deg2arcsec/RAscale
    RAmax = args.x_center + 0.5*args.width*args.pixel_scale/deg2arcsec/RAscale
    DECmin = args.y_center - 0.5*args.height*args.pixel_scale/deg2arcsec
    DECmax = args.y_center + 0.5*args.height*args.pixel_scale/deg2arcsec
    
    # Calculate margin size in degrees (sources outside of our image margins
    # are always skipped, for speed, even if their tails might overlap our image)
    margin = args.margin/deg2arcmin
    
    # Convert the band into an integer index
    bandIndex = "ugrizy".find(args.band)

    # Calculate the sky noise level in stacked ADU / pixel
    skyNoise = math.sqrt(2*args.nvisits*args.sky_level)
    
    # Calculate the stacked pixel ADU threshold cut to use
    pixelCut = args.sn_cut*skyNoise

    # Calculate the corresponding surface brightness cut to use
    sbCut = pixelCut/(args.pixel_scale*args.pixel_scale)
    
    print 'Simulating %s-band observations with flux(AB24) = %.3f ADU/pixel/exp for typical galaxy.' %(
        args.band,args.flux_norm)
    print 'Simulating %d visits with stacked sky noise level %.3f ADU/pixel.' % (
        args.nvisits,skyNoise)
    print 'Will keep all stacked pixels > %.3f ADU (%.1f ADU/arcsec^2)' % (pixelCut,sbCut)

    # Initialize finite difference calculations if necessary
    if args.partials:
        args.stamps = True
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

    # Open the source input catalog to use and read the header line
    cat = open(args.input)
    cathdr = cat.readline().split()
    logger.info('Reading input catalog %r with fields:\n%s' % (args.input,','.join(cathdr)))

    # Open the output catalog to write
    outname = args.output + '_catalog.dat'
    outcat = open(outname,'w')
    if args.verbose:
        logger.info('Creating output catalog to %r' % outname)
    
    # Initialize the list of per-object stamp HDUs we will fill
    hdu = pyfits.PrimaryHDU()
    hduList = pyfits.HDUList([hdu])

    # Loop over catalog entries
    nkeep = lineno = 0
    for line in cat:
        lineno += 1

        if args.only_line > 0 and lineno != args.only_line:
            continue

        # position on the sky in degrees
        cols = line.split()
        RA = float(cols[1])
        DEC = float(cols[2])
        
        # skip sources outside our margins
        if RA < RAmin-margin or RA > RAmax+margin or DEC < DECmin-margin or DEC > DECmax+margin:
            continue
        
        # Look up source AB magnitude in the requested band
        abMag = float(cols[19+bandIndex])
        # Calculate total flux in ADU units
        flux = args.flux_norm*math.pow(10,24-abMag)
        # Scale flux to number of vists (extra factor of 2 because 1 visit = 2 exposures)
        flux = 2*args.nvisits*flux
        # Skip objects whose total flux is below our pixel threshold
        if flux < pixelCut:
            continue
        
        # Look up the disk and bulge fluxes, which are provided in the catalog as
        # color-independent magnitudes.
        bulgeMag = float(cols[9])
        diskMag = float(cols[10])
        if bulgeMag > 0 and diskMag > 0:
            bulgeFraction = 1./(1.+math.exp(-(diskMag-bulgeMag)/2.5))
        elif bulgeMag > 0:
            bulgeFraction = 1
        else:
            bulgeFraction = 0
        if args.no_bulge:
            if bulgeFraction == 1:
                # skip sources that are only bulge
                continue
            else:
                bulgeFraction = 0
        
        # Calculate bounding-box padding in arcsecs for pixel and psf convolution
        (w_pad,h_pad) = moffatBounds(args.psf_beta,flux,args.psf_fwhm,1,0,sbCut)
        w_pad += args.pixel_scale
        h_pad += args.pixel_scale
        
        # Get disk component parameters
        if bulgeFraction < 1:
            hlr_d = float(cols[6]) # in arcsecs
            if hlr_d <= 0:
                raise RuntimeError('Unexpected hlr_d <= 0')
            pa_d = float(cols[8]) # position angle in degrees
            a_d = float(cols[16]) # major axis length in arcsecs
            b_d = float(cols[18]) # minor axis length in arcsecs
            # Calculate sheared ellipse aspect ratio
            q_d = b_d/a_d # between 0.2 and 1
            # Convert position angle from degrees to radians
            pa_d = pa_d*deg2rad
            # Calculate bounding box in arcsecs without psf or pixel convolution
            (w_d,h_d) = sersicBounds(1,flux,hlr_d,q_d,pa_d,sbCut)
            # Add padding for psf and pixel convolution, unless disk's surface
            # density is always below (1-bulgeFraction)*sbCut
            if (w_d,h_d) != (0,0):
                w_d += w_pad
                h_d += h_pad
        else:
            (hlr_d,q_d,pa_d) = (0,0,0)
            (w_d,h_d) = (0,0)
        
        # Get bulge component parameters
        if bulgeFraction > 0:
            hlr_b = float(cols[5]) # in arcsecs
            if hlr_b <= 0:
                raise RuntimeError('Unexpected hlr_b <= 0')
            pa_b = float(cols[7]) # position angle in degrees
            a_b = float(cols[15]) # major axis length in arcsecs
            b_b = float(cols[17]) # minor axis length in arcsecs
            # Calculate sheared ellipse aspect ratio
            q_b = b_b/a_b # between 0.2 and 1
            # Convert position angle from degrees to radians
            pa_b = pa_b*deg2rad
            # Calculate bounding box in arcsecs without psf or pixel convolution
            (w_b,h_b) = sersicBounds(4,flux,hlr_b,q_b,pa_b,sbCut)
            # Add padding for psf and pixel convolution, unless bulge's surface
            # density is always below bulgeFraction*sbCut
            if (w_b,h_b) != (0,0):
                w_b += w_pad
                h_b += h_pad
        else:
            (hlr_b,q_b,pa_b) = (0,0,0)
            (w_b,h_b) = (0,0)
        
        # Combined the bulge and disk bounding boxes
        width = max(w_d,w_b)
        height = max(h_d,h_b)
        
        # Skip this source if its pixels would all be below pixelCut
        if (width,height) == (0,0):
            continue
        
        # Calculate the offsets of this source from our image's bottom left corner in pixels
        # (which might be negative, or byeond our image bounds)
        xoffset = (RA - RAmin)*deg2arcsec/args.pixel_scale*RAscale
        yoffset = (DEC - DECmin)*deg2arcsec/args.pixel_scale
        
        # Calculate the integer coordinates of the image pixel that contains the source center
        # (using the convention that the bottom left corner pixel has coordinates 1,1)
        xpixels = int(math.ceil(xoffset))
        ypixels = int(math.ceil(yoffset))

        # Calculate the stamp size to use as width = 2*xhalf+1 and height = 2*yhalf+1.
        # We always round up to an odd integer so that flux is consistently centered
        # (see Issue #380).
        xhalf = int(math.ceil(width/args.pixel_scale))
        yhalf = int(math.ceil(height/args.pixel_scale))
        
        # Trim the stamp so that the source is still centered but we do not extend
        # beyond the final field image. This will only trim pixels above pixelCut
        # that lie outside the field.
        if xpixels-xhalf < 1 and xpixels+xhalf > args.width:
            xhalf = max(xpixels-1,args.width-xpixels)
        if ypixels-yhalf < 1 and ypixels+yhalf > args.height:
            yhalf = max(ypixels-1,args.height-ypixels)

        # Build this source's stamp bounding box
        bbox = galsim.BoundsI(xpixels-xhalf,xpixels+xhalf,ypixels-yhalf,ypixels+yhalf)
        
        # Skip objects that don't overlap our field
        if (bbox & field.bounds).area() == 0:
            continue

        # If we get this far, we are definitely rendering this source (but it might
        # still get trimmed out later)
        nkeep += 1
        logger.info('Rendering stamp %d (line %d) with w x h = %d x %d' %
            (nkeep,lineno,2*xhalf+1,2*yhalf+1))

        # Calculate the pixel coordinates of the stamp center.
        xstamp = 0.5*(bbox.xmin + bbox.xmax)
        ystamp = 0.5*(bbox.ymin + bbox.ymax)

        # Calculate the subpixel shift in arcsecs (not pixels!) of the source center
        # relative to the stamp center. Note that the resulting shift may be more than
        # one pixel in either direction because of the clipping operation above.
        xshift = (xoffset - (xstamp-0.5))*args.pixel_scale
        yshift = (yoffset - (ystamp-0.5))*args.pixel_scale

        if args.verbose:
            logger.info('    flux: %.3g ADU (%s-band AB %.1f)' % (flux,args.band,abMag))
            logger.info('  bounds: [%d:%d,%d:%d] pixels' % (bbox.xmin,bbox.xmax,bbox.ymin,bbox.ymax))
            logger.info('   shift: (%f,%f) arcsecs = (%f,%f) pixels' %
                (xshift,yshift,xshift/args.pixel_scale,yshift/args.pixel_scale))
            logger.info('    disk: f = %f, hlr = %f arcsec, q = %f, beta = %f rad' %
                (1-bulgeFraction,hlr_d,q_d,pa_d))
            logger.info('   bulge: f = %f, hlr = %f arcsec, q = %f, beta = %f rad' %
                (bulgeFraction,hlr_b,q_b,pa_b))
            logger.info('    bbox: disk (%.1f,%.1f) bulge (%.1f,%.1f) pad (%.1f,%.1f) pixels' %
                (w_d,h_d,w_b,h_b,w_pad,h_pad))
        
        # Define the nominal source parameters for rendering this object within its stamp
        params = {
            'flux':flux, 'bulgeFraction': bulgeFraction,
            'xc':xshift, 'yc':yshift,
            'hlr_d':hlr_d, 'q_d':q_d, 'beta_d':pa_d,
            'hlr_b':hlr_b, 'q_b':q_b, 'beta_b':pa_b,
            'g1':args.g1, 'g2': args.g2
        }

        # Create stamps for the galaxy with and w/o the psf applied
        gal = createSource(**params)
        nopsf = createStamp(gal,None,pix,bbox)
        nominal = createStamp(gal,psf,pix,bbox)        
        mask = createMask(nominal,pixelCut,args)
        trimmed = mask.bounds
        if not args.no_trim and args.verbose:
            logger.info(' trimmed: [%d:%d,%d:%d] pixels' % (trimmed.xmin,trimmed.xmax,trimmed.ymin,trimmed.ymax))
        
        # Add the nominal galaxy to the full field image
        overlap = nominal.bounds & field.bounds
        field[overlap] += nominal[overlap]
        
        # Initialize the datacube of stamps that we will save for this object
        datacube = [ ]
        singleExposureNorm=1./(2*args.nvisits)
        saveStamp(datacube,singleExposureNorm*nopsf,trimmed,args)
        saveStamp(datacube,singleExposureNorm*nominal,trimmed,args)
        if not saveStamp(datacube,mask,trimmed,args):
            # this stamp's mask falls completely outside our field
            logger.info('*** stamp %d not saved' % nkeep)
            nkeep -= 1
            continue

        if args.partials:
            # Specify the amount to vary each parameter for partial derivatives
            # (we don't use a dictionary here since we want to control the order)
            variations = [
                ('xc',args.pixel_scale/3.), ('yc',args.pixel_scale/3.),
                ('hlr_d',0.05*hlr_d),('hlr_b',0.05*hlr_b),
                ('g1',0.03), ('g2',0.03)
            ]
            # loop over source parameters to vary
            for (pname,delta) in variations:
                # create stamps for each variation of this parameter
                newparams = params.copy()
                partial = galsim.ImageD(bbox)
                partial.setScale(pix.getXWidth())
                # delta might be zero, e.g., for hlr_b when bulge fraction = 0
                if delta > 0:
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
                assert saveStamp(datacube,singleExposureNorm*partial,trimmed,args)

        # Add a new HDU with a datacube for this object's stamps
        # We don't use compression = 'gzip_tile' for now since it is lossy
        # and mathematica cannot Import it.
        galsim.fits.writeCube(datacube, hdu_list = hduList)

        # Write an entry for this object to the output catalog
        print >>outcat, lineno,xstamp,ystamp,abMag,flux/(2*args.nvisits)

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
    if args.stamps:
        outname = args.output + '_stamps.fits'
        logger.info('Saving %d stamps to %r' % (nkeep,outname))
        galsim.fits.writeFile(outname, hduList)

    # Close catalog files
    cat.close()
    outcat.close()

if __name__ == "__main__":
    main()
