#!/usr/bin/env python
#######################################################################################
## Created by David Kirkby, University of California, Irvine <dkirkby@uci.edu>
#######################################################################################

import sys
import os
import math
import argparse

import logging
import galsim
import pyfits
import numpy

twopi = 2*math.pi
deg2rad = math.pi/180.
deg2arcsec = 3600.
deg2arcmin = 60.

"""
Returns a (disk,bulge) tuple of source objects using the specified parameters.
"""
def createSource(flux,bulgeFraction,xc,yc,hlr_d,q_d,beta_d,hlr_b,q_b,beta_b,g1,g2):
    # Define the disk component, if any
    if bulgeFraction < 1:
        disk = galsim.Exponential(flux = flux*(1-bulgeFraction), half_light_radius = hlr_d)
        # Apply intrinsic shear
        disk.applyShear(q = q_d, beta = beta_d*galsim.radians)
        # Apply cosmic shear
        disk.applyShear(g1 = g1, g2 = g2)
        # Shift to this object's centroid
        disk.applyShift(dx = xc, dy = yc)
    else:
        disk = None
    # Define the bulge component, if any
    if bulgeFraction > 0:
        bulge = galsim.DeVaucouleurs(flux = flux*bulgeFraction, half_light_radius = hlr_b)
        # Apply intrinsic shear
        bulge.applyShear(q = q_b, beta = beta_b*galsim.radians)
        # Apply cosmic shear
        bulge.applyShear(g1 = g1, g2 = g2)
        # Shift to this object's centroid
        bulge.applyShift(dx = xc, dy = yc)
    else:
        bulge = None
    return (disk,bulge)

"""
Renders the specified source convolved with a psf (which might be None)
and pixel response into a postage stamp with the specified bounding box.
"""
def renderStamp(src,psf,pix,bbox):
    stamp = galsim.ImageD(bbox)
    stamp.setScale(pix.getXWidth())
    if src:
        if psf == None:
            obj = galsim.Convolve([src,pix], real_space = True)
        else:
            gsp = galsim.GSParams(maximum_fft_size=16384)
            obj = galsim.Convolve([src,psf,pix],gsparams=gsp)
        obj.draw(image = stamp)
    return stamp

def createStamp(src,psf,pix,bbox):
    (disk,bulge) = src
    diskStamp = renderStamp(disk,psf,pix,bbox)
    bulgeStamp = renderStamp(bulge,psf,pix,bbox)
    return diskStamp + bulgeStamp

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
0 < q <= 1 and moffatBeta > 1 are dimensionless. The returned (dx,dy) are in
arcsecs. See boundingBox above for details.
"""
def moffatBounds(moffatBeta,flux,fwhm,q,beta,f0):
    # Check that moffatBeta is valid
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
Note that if all pixels are below threshold, then the returned mask will
contain only the central pixel with image.array.sum() == 0.
"""
def createMask(image,threshold,args):
    # create an empty mask image with the same dimensions as the input image
    box = image.bounds
    mask = galsim.ImageS(box)
    mask.setScale(image.getScale())
    borderMax = 0.
    lastRow = box.ymax - box.ymin
    lastPixel = box.xmax - box.xmin
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
            if rowIndex == 0 or rowIndex == lastRow or pixelIndex == 0 or pixelIndex == lastPixel:
                # update the largest pixel value on our 1-pixel wide border
                borderMax = max(borderMax,pixelValue)
            if pixelValue >= threshold:
                mask.array[rowIndex,pixelIndex] = 1
                if not args.no_trim:
                    xmin = min(x,xmin)
                    xmax = max(x,xmax)
                    ymin = min(y,ymin)
                    ymax = max(y,ymax)
    # is the stamp too small to contain the threshold contour?
    if borderMax > threshold:
        print '### stamp truncated at %.1f > %.1f ADU' % (borderMax,threshold)
        # build a new mask using the border max as the threshold
        return createMask(image,borderMax,args)
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

def getPsfBoundsEstimator(psf,pix,size):
    stamp = galsim.ImageD(2*size,2*size)
    stamp.setScale(pix.getXWidth())
    obj = galsim.Convolve([psf,pix])
    obj.draw(image=stamp)
    # build the circularized psf profile
    profile = numpy.zeros(size,dtype=float)
    for x in range(2*size):
        for y in range(2*size):
            dx = x - size + 0.5
            dy = y - size + 0.5
            r = math.sqrt(dx*dx + dy*dy)
            ipix = int(math.floor(r))
            if ipix < size:
                profile[ipix] = max(profile[ipix],stamp.array[x,y])
    # Create a function that gives the size of bounding box necessary to contain
    # psf pixels down to the specified threshold assuming the specified total flux.
    # The return value is clipped at 2*size for large fluxes.
    def estimator(flux,threshold):
        index = 0
        while index < size:
            if flux*profile[index] < threshold:
                return 2*index+1
            index += 1
        return 2*size
    return estimator

# Returns the combined size and ellipticity for the specified disk and bulge components,
# assuming they have the same centroid.
def combineEllipticities(hlr_d,q_d,pa_d,hlr_b,q_b,pa_b,f_b):
    # ensure that single-component models give correct results
    if f_b == 0:
        q_b = 1
    elif f_b == 1:
        q_d = 1

    # calculate the disk and bulge component ellipticities
    ed = (1-q_d)/(1+q_d)
    ed1 = ed*math.cos(2*pa_d)
    ed2 = ed*math.sin(2*pa_d)
    eb = (1-q_b)/(1+q_b)
    eb1 = eb*math.cos(2*pa_b)
    eb2 = eb*math.sin(2*pa_b)

    # calculate the corresponding second-moment tensors assuming unit total flux
    cd = 2.13003
    nd = cd*(hlr_d/(1-ed*ed))**2
    Qd11 = nd*(1+ed*ed+2*ed1)
    Qd12 = nd*2*ed2
    Qd22 = nd*(1+ed*ed-2*ed1)
    detQd = Qd11*Qd22 - Qd12*Qd12
    cb = 21.6793
    nb = cb*(hlr_b/(1-eb*eb))**2
    Qb11 = nb*(1+eb*eb+2*eb1)
    Qb12 = nb*2*eb2
    Qb22 = nb*(1+eb*eb-2*eb1)
    detQb = Qb11*Qb22 - Qb12*Qb12

    # add the component second-moment tensors
    Q11 = (1-f_b)*Qd11 + f_b*Qb11
    Q12 = (1-f_b)*Qd12 + f_b*Qb12
    Q22 = (1-f_b)*Qd22 + f_b*Qb22
    detQ = Q11*Q22 - Q12*Q12
    size = math.pow(detQ,0.25)

    # calculate the corresponding combined ellipticity
    denom = Q11 + Q22 + 2*math.sqrt(detQ)
    e1 = (Q11 - Q22)/denom
    e2 = 2*Q12/denom

    """
    # check direct calculation of emag when pa_d == pa_b
    emag = math.sqrt(e1*e1 + e2*e2)
    wd = (1-f_b)*cd*hlr_d**2*(1+q_d)**2/(8*q_d**2) if f_b < 1 else 0
    wm = f_b*cb*hlr_b**2*(1+q_b)**2/(8*q_b**2) if f_b > 0 else 0
    ep = wd*(1+q_d**2) + wm*(1+q_b**2)
    em = wd*(1-q_d**2) + wm*(1-q_b**2)
    emag2 = em/(ep+math.sqrt(ep*ep-em*em))
    print 'emag:',emag-emag2
    """

    return (size,e1,e2)

def main():

    # Parse command-line args
    parser = argparse.ArgumentParser()
    parser.add_argument("-v", "--verbose", action = "store_true",
        help = "provide more verbose output")
    parser.add_argument("-i", "--input", default = 'gcat.dat',
        help = "name of input catalog to read")
    parser.add_argument("-o","--output", default = 'catout',
        help = "base name of output files to write")
    parser.add_argument("--catscan", action = "store_true",
        help = "build output catalog only, with no rendering")
    parser.add_argument("-x","--x-center", type = float, default = 0.5,
        help = "central RA of image (degrees)")
    parser.add_argument("-y","--y-center", type = float, default = 0.0,
        help = "central DEC of image (degrees)")
    parser.add_argument("--width", type = int, default = 512,
        help = "image width (pixels)")
    parser.add_argument("--height", type = int, default = 512,
        help = "image height (pixels)")
    parser.add_argument("--max-size", type = float, default = 20.,
        help = "flux from any object is truncated beyond this size (arcsecs)")
    parser.add_argument("--no-margin", action = "store_true",
        help = "do not simulate the tails of objects just outside the field")
    parser.add_argument("--pixel-scale", type = float, default = 0.2,
        help = "pixel scale (arscecs/pixel)")
    parser.add_argument("--psf-fwhm", type = float, default = 0.7,
        help = "psf full-width-half-max in arcsecs (zero for no psf)")
    parser.add_argument("--psf-beta", type = float, default = 0.0,
        help = "psf Moffat parameter beta (uses Kolmogorov psf if beta <= 0)")
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
    parser.add_argument("--no-disk", action = "store_true",
        help = "do not include any galactic disk (Sersic n=1) components")
    parser.add_argument("--no-bulge", action = "store_true",
        help = "do not include any galactic bulge (Sersic n=4) components")
    parser.add_argument("--shape", action = "store_true",
        help = "run HSM adaptive moments calculation on no-psf stamp")
    parser.add_argument("--render-nopsf", action = "store_true",
        help = "save a stamp rendered without any psf for each galaxy")
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
        if args.psf_beta > 0:
            psf = galsim.Moffat(beta = args.psf_beta, fwhm = args.psf_fwhm)
        else:
            psf = galsim.Kolmogorov(fwhm = args.psf_fwhm)
        psfBounds = getPsfBoundsEstimator(psf,pix,int(math.ceil(0.5*args.max_size/args.pixel_scale)))
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
    if args.no_margin:
        margin = 0
    else:
        margin = 0.5*args.max_size/deg2arcsec
    
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
    print 'Simulating %d visits with stacked sky noise level %.3f ADU/pixel (%.3f sky ADU/pixel/exp)' % (
        args.nvisits,skyNoise,args.sky_level)
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
    catoutname = args.output + '_catalog.dat'
    outcat = open(catoutname,'w')
    
    # Initialize the list of per-object stamp HDUs we will fill
    hdu = pyfits.PrimaryHDU()
    hduList = pyfits.HDUList([hdu])

    # Loop over catalog entries
    ncat = nkeep = lineno = 0
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
        
        # Calculate the offsets of this source from our image's bottom left corner in pixels
        # (which might be negative or byeond our image bounds because of the margins)
        xoffset = (RA - RAmin)*deg2arcsec/args.pixel_scale*RAscale
        yoffset = (DEC - DECmin)*deg2arcsec/args.pixel_scale
        
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
            bulgeFraction = 1./(1.+math.pow(10,-(diskMag-bulgeMag)/2.5))
        elif bulgeMag > 0:
            bulgeFraction = 1
        else:
            bulgeFraction = 0
        if args.no_bulge:
            if bulgeFraction == 1:
                # skip sources that are only bulge
                continue
            else:
                # only render flux associated with the disk
                flux = (1-bulgeFraction)*flux
                bulgeFraction = 0
        elif args.no_disk:
            if bulgeFraction == 0:
                # skip sources that are only disk
                continue
            else:
                # only render flux associated with the bulge
                flux = bulgeFraction*flux
                bulgeFraction = 1
        
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
            if not args.catscan:
                (w_d,h_d) = sersicBounds(1,flux,hlr_d,q_d,pa_d,sbCut)
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
            if not args.catscan:
                (w_b,h_b) = sersicBounds(4,flux,hlr_b,q_b,pa_b,sbCut)
        else:
            (hlr_b,q_b,pa_b) = (0,0,0)
            (w_b,h_b) = (0,0)

        # Combine the bulge and disk ellipticities
        (size,e1,e2) = combineEllipticities(hlr_d,q_d,pa_d,hlr_b,q_b,pa_b,bulgeFraction)

        # Write an entry for this object to the output catalog
        print >>outcat, lineno,xoffset,yoffset,abMag,flux/(2*args.nvisits),size,e1,e2
        ncat += 1

        # All done now in catscan mode
        if args.catscan:
            continue

        # Combine the bulge and disk bounding boxes
        width = max(w_d,w_b)
        height = max(h_d,h_b)

        # Estimate the (round) bounding box for the psf in arscecs given our total flux
        psfSize = psfBounds(flux,pixelCut)*args.pixel_scale if psf else 0
        # Add the psf size in quadrature
        width = math.sqrt(width*width + psfSize*psfSize)
        height = math.sqrt(height*height + psfSize*psfSize)

        # Truncate the bounding box, if necessary
        if width > args.max_size or height > args.max_size:
            logger.info('...truncating bbbox from (%.1f,%.1f)' % (width,height))
            width = min(width,args.max_size)
            height = min(height,args.max_size)
        
        # Skip this source if its pixels would all be below pixelCut (can this ever happen?)
        if (width,height) == (0,0):
            continue
        
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
            logger.info('    bbox: disk (%.1f,%.1f) bulge (%.1f,%.1f) psf %.1f arcsec' %
                (w_d,h_d,w_b,h_b,psfSize))
            logger.info('    size: %.2f pixels' % size)
            logger.info('   shear: (g1,g2) = (%.6f,%.6f)' % (e1,e2))
        
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
        if args.render_nopsf:
            # We don't do this by default for speed
            nopsf = createStamp(gal,None,pix,bbox)
        else:
            # Use an empty placeholder so we don't change the shape of the output
            nopsf = galsim.ImageD(bbox)
        nominal = createStamp(gal,psf,pix,bbox)

        # Create a mask for pixels above threshold
        mask = createMask(nominal,pixelCut,args)
        if mask.array.sum() == 0:
            # this stamp has no pixels above threshold
            logger.info('*** stamp %d is below threshold' % nkeep)
            nkeep -= 1
            continue
        trimmed = mask.bounds
        if not args.no_trim and args.verbose:
            logger.info(' trimmed: [%d:%d,%d:%d] pixels' % (trimmed.xmin,trimmed.xmax,trimmed.ymin,trimmed.ymax))

        # Add the nominal galaxy to the full field image after applying the threshold mask
        # (the mask must be the second term in the product so that the result is double precision)
        masked = nominal[trimmed]*mask
        overlap = trimmed & field.bounds
        if overlap.area() == 0:
             # this stamp's mask falls completely outside our field
            logger.info('*** stamp %d does not overlap field' % nkeep)
            nkeep -= 1
            continue
        field[overlap] += masked[overlap]
        
        # Initialize the datacube of stamps that we will save for this object
        datacube = [ ]
        singleExposureNorm=1./(2*args.nvisits)
        assert saveStamp(datacube,singleExposureNorm*nopsf,trimmed,args)
        assert saveStamp(datacube,singleExposureNorm*nominal,trimmed,args)
        assert saveStamp(datacube,mask,trimmed,args)

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

    # Write the full field image to a separate file
    if not args.catscan:
        outname = args.output + '_field.fits'
        logger.info('Saving full field to %r' % outname)
        galsim.fits.write(field,outname)

    # Write out the full field image with noise added
    if not args.catscan and args.sky_level > 0:
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
    logger.info('Wrote %d of %d catalog entries to %r' % (ncat,lineno,catoutname))

if __name__ == "__main__":
    main()
