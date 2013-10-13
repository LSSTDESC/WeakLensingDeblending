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

def createComponent(type,electrons,xc,yc,hlr,q,beta,g1,g2):
    # create a radial profile of the requested type and size
    comp = type(flux = electrons, half_light_radius = hlr)
    # set the intrinsic shape
    comp.applyShear(q = q, beta = beta*galsim.radians)
    # apply cosmic shear
    comp.applyShear(g1 = g1, g2 = g2)
    # shift to this component's centroid
    comp.applyShift(dx = xc, dy = yc)
    return comp

"""
Returns a (disk,bulge) tuple of source objects using the specified parameters.
Note that f_d and f_b are fractions of the total flux and need not sum to one.
"""
def createSource(
    total_flux,f_d,f_b,
    x_d,y_d,hlr_d,q_d,beta_d,
    x_b,y_b,hlr_b,q_b,beta_b,
    dx,dy,relsize,dbeta,
    g1,g2):
    # Define the disk component, if any
    if f_d > 0:
        disk = createComponent(galsim.Exponential,
            total_flux*f_d,x_d+dx,y_d+dy,hlr_d*relsize,q_d,beta_d+dbeta,g1,g2)
    else:
        disk = None
    # Define the bulge component, if any
    if f_b > 0:
        bulge = createComponent(galsim.DeVaucouleurs,
            total_flux*f_b,x_b+dx,y_b+dy,hlr_b*relsize,q_b,beta_b+dbeta,g1,g2)
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
The input flux should be in electrons, hlr in arscecs, beta in radians, f0 in
elec/arcsec^2. 0 < q <= 1 is dimensionless. The returned (dx,dy) are in
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
should be in electrons, fwhm in arcsecs, beta in radians, f0 in elec/arcsec^2.
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
input image pixel value is above or below the specified threshold in electrons.
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
                xmin = min(x,xmin)
                xmax = max(x,xmax)
                ymin = min(y,ymin)
                ymax = max(y,ymax)
    # is the stamp too small to contain the threshold contour?
    if borderMax > threshold:
        print '### stamp truncated at %.1f > %.1f electrons' % (borderMax,threshold)
        # build a new mask using the border max as the threshold
        return createMask(image,borderMax,args)
    trimmed = galsim.BoundsI(xmin,xmax,ymin,ymax)
    mask = mask[trimmed]
    return mask

"""
Performs any final processing on stamp, controlled by args, then appends it to stamps.
Returns True if the stamp was saved, or otherwise False.
"""
def saveStamp(stamps,stamp,args):
    # Clip the stamp so that does not extend beyond the field image. This results
    # in potentially smaller files with sources that might not be centered.
    if not args.no_clip:
        overlap = stamp.bounds & galsim.BoundsI(1,args.width,1,args.height)
        if overlap.area() == 0:
            # skip this stamp if it falls completely outside our field (after trimming)
            return False
        stamp = stamp[overlap]
    # Convert normalization from elec/exposure to elec/second
    stamp = stamp/args.exposure_time
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
    cd = 1.06502
    nd = cd*(hlr_d/(1-ed*ed))**2
    Qd11 = nd*(1+ed*ed+2*ed1)
    Qd12 = nd*2*ed2
    Qd22 = nd*(1+ed*ed-2*ed1)
    detQd = Qd11*Qd22 - Qd12*Qd12
    cb = 10.8396
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
    #semiMajorAxis = math.sqrt(0.5*(Q11+Q22+math.sqrt((Q11-Q22)**2+4*Q12**2)))

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

def signalToNoiseRatio(stamp,pixelNoise):
    flat = stamp.array.reshape(-1)
    snr = math.sqrt(numpy.dot(flat,flat)/pixelNoise)
    return snr

# Returns True if the stamps s1 and s2 have overlapping pixels with non-zero flux.
def overlapping(s1,s2):
    # test for overlapping bounding boxes
    overlapBounds = s1.bounds & s2.bounds
    if overlapBounds.area() == 0:
        return False
    # test for overlapping flux within the overlapping pixels
    overlapFluxProduct = numpy.sum(s1[overlapBounds].array * s2[overlapBounds].array)
    return False if overlapFluxProduct == 0 else True

# Assigns a group ID to each stamp in stamps based on its overlaps with other stamps.
def analyzeOverlaps(stamps):
    groupID = range(len(stamps))
    groupSize = [1]*len(stamps)
    for (i1,s1) in enumerate(stamps):
        for (i2,s2) in enumerate(stamps[:i1]):
            if overlapping(s1,s2):
                # get the current group IDs of these overlapping stamps
                gid1 = groupID[i1]
                gid2 = groupID[i2]
                if gid1 == gid2:
                    continue
                # decide which group joins the other
                gnew = min(gid1,gid2)
                gold = max(gid1,gid2)
                # re-assign all stamps in gold to gnew
                for i in range(i1+1):
                    if groupID[i] == gold:
                        groupID[i] = gnew
                        groupSize[gnew] += 1
                        groupSize[gold] -= 1
    return (groupID,groupSize)

# Builds the Fisher matrix from the specified array of npar*(npar+1)/2 Fisher images and
# calculates the corresponding shape-measurment error, if possible.
def shapeError(npar,fisherImages,mask):
    # calculate the Fisher matrix elements by summing pixels of the specified Fisher matrix images
    fisherMatrix = numpy.zeros((npar,npar))
    index = 0
    for i in range(npar):
        for j in range(i,npar):
            fisherMatrix[i,j] = numpy.sum(fisherImages[index]*mask)
            if i != j:
                fisherMatrix[j,i] = fisherMatrix[i,j]
            index += 1
    # try to calculate corresponding shape measurement error, which will fail unless
    # the Fisher matrix is invertible
    try:
        fullCov = numpy.linalg.inv(fisherMatrix)
        # this is where we assume that the last 2 variations are g1,g2
        varEps = 0.5*(fullCov[-2,-2]+fullCov[-1,-1])
        # variance might be negative if inverse has large numerical errors
        sigmaEps = 0 if varEps <= 0 else math.sqrt(varEps)
    except numpy.linalg.LinAlgError:
        # assign a shape-measurement error of zero if the Fisher matrix is not invertible.
        sigmaEps = 0.
    return sigmaEps

# Calculate shape measurment errors with the specified purity cuts. Returns a tuple of the
# corresponding errors, in a list, and an integer-valued image that identifies the purity
# regions by assigning each pixel the value of the largest index such that
# nominal > purity[index]*field (or zero if this criteria is not met for any purity).
def shapeErrorsAnalysis(npar,nominal,fisherImages,field,purities):
    # find the overlap of this object in the full field
    overlap = nominal.bounds & field.bounds
    subNominal = nominal[overlap].array
    subField = field[overlap].array
    regions = galsim.ImageI(nominal.bounds)
    regionsArray = regions.array
    errors = [ ]
    for (i,purity) in enumerate(purities):
        mask = (subNominal > purity*subField)
        regionsArray = numpy.maximum(regionsArray,i*mask)
        sigeps = shapeError(npar,fisherImages,mask)
        errors.append(sigeps)
    regions.array[:] = regionsArray[:]
    return (errors,regions)

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
    parser.add_argument("--max-size", type = float, default = 20.,
        help = "flux from any object is truncated beyond this size (arcsecs)")
    parser.add_argument("--no-margin", action = "store_true",
        help = "do not simulate the tails of objects just outside the field")
    parser.add_argument("--pixel-scale", type = float, default = 0.2,
        help = "pixel scale (arscecs/pixel)")
    parser.add_argument("--airmass", type = float, default = 1.2,
        help = "airmass value to use for atmospheric PSF and extinction")
    parser.add_argument("--extinction", type = float, default = 0.07,
        help = "atmospheric extinction coefficient")
    parser.add_argument("--zenith-fwhm", type = float, default = 0.67,
        help = "atmospheric psf full-width-half-max in arcsecs at zenith")
    parser.add_argument("--instrumental-fwhm", type = float, default = 0.4,
        help = "instrumental psf full-width-half-max in arcsecs")
    parser.add_argument("--psf-beta", type = float, default = 0.0,
        help = "psf Moffat parameter beta (uses Kolmogorov psf if beta <= 0)")
    parser.add_argument("--band", choices = ['u','g','r','i','z','y'], default = 'i',
        help = "LSST imaging band to use for source fluxes")
    parser.add_argument("--zero-point", type = float, default = 41.5,
        help = "zero point for converting magnitude to detected signal in elec/sec")
    parser.add_argument("--sky-brightness", type = float, default = 20.0,
        help = "sky brightness in mag/sq.arcsec.")
    parser.add_argument("--sn-cut", type = float, default = 0.5,
        help = "keep all pixels above this signal-to-noise ratio cut")
    parser.add_argument("--exposure-time", type = float, default = 6900.,
        help = "full-depth exposure time in seconds")
    parser.add_argument("--g1", type = float, default = 0.,
        help = "constant shear component g1 to apply")
    parser.add_argument("--g2", type = float, default = 0.,
        help = "constant shear component g2 to apply")
    parser.add_argument("--save-field", action = "store_true",
        help = "save full field image without noise")
    parser.add_argument("--save-noise", action = "store_true",
        help = "save full field image with random noise added")
    parser.add_argument("--stamps", action = "store_true",
        help = "save postage stamps for each source (normalized to 1 exposure)")
    parser.add_argument("--no-clip", action = "store_true",
        help = "do not clip stamps to the image bounds")
    parser.add_argument("--no-disk", action = "store_true",
        help = "do not include any galactic disk (Sersic n=1) components")
    parser.add_argument("--no-bulge", action = "store_true",
        help = "do not include any galactic bulge (Sersic n=4) components")
    parser.add_argument("--shape", action = "store_true",
        help = "run HSM adaptive moments calculation on no-psf stamp")
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
    atmos_fwhm = args.zenith_fwhm*math.pow(args.airmass,0.6)
    fwhm = math.sqrt(atmos_fwhm**2 + args.instrumental_fwhm**2)
    logger.info('Using PSF fwhm = %.4f" (%.4f" zenith => %.4f" at X = %.3d, %.4f" instrumental)' %
        (fwhm,args.zenith_fwhm,atmos_fwhm,args.airmass,args.instrumental_fwhm))
    if fwhm > 0:
        if args.psf_beta > 0:
            psf = galsim.Moffat(beta = args.psf_beta, fwhm = fwhm)
        else:
            psf = galsim.Kolmogorov(fwhm = fwhm)
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
    
    # Calculate the sky background rate in elec/sec/pixel
    skyRate = args.zero_point*math.pow(10,-0.4*(args.sky_brightness-24))*args.pixel_scale**2

    # Calculate the mean sky noise level for the full exposure time in elec/pixel
    skyNoise = math.sqrt(args.exposure_time*skyRate)
    
    # Calculate the pixel threshold cut to use in detected electrons during the full exposure
    pixelCut = args.sn_cut*skyNoise

    # Calculate the corresponding surface brightness cut to use
    sbCut = pixelCut/(args.pixel_scale*args.pixel_scale)
    
    print 'Simulating %s-band observations with AB24 zero point %.3f elec/sec, sky rate = %.3f elec/sec/pixel' %(
        args.band,args.zero_point,skyRate)
    print 'Simulating %.1fs exposure with total sky noise level %.3f elec/pixel (%.3f mag/sq.arcsec.)' % (
        args.exposure_time,skyNoise,args.sky_brightness)
    print 'Will keep all stacked pixels > %.3f elec (%.1f elec/arcsec^2)' % (pixelCut,sbCut)

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

    # Open the source input catalog to use and initialize a keyword-based lookup for catalog entries
    cat = open(args.input)
    catFields = cat.readline().split()
    catDict = dict(zip(catFields,range(len(catFields))))
    if args.verbose:
        logger.info('Reading input catalog %r with fields:\n%s' % (args.input,','.join(catFields)))

    # Initialize the output catalog in memory
    outputCatalog = [ ]
    
    # Initialize the list of per-object stamp HDUs we will fill
    hdu = pyfits.PrimaryHDU()
    hduList = pyfits.HDUList([hdu])
    stampList = [ ]
    fisherImagesList = [ ]
    nvar = 0 # declared here so it stays in scope after loop over galaxies

    # Loop over catalog entries
    nkeep = lineno = 0
    for line in cat:
        lineno += 1

        if args.only_line > 0 and lineno != args.only_line:
            continue

        # prepare to read this catalog entry
        entryCols = line.split()
        def catalog(fieldName,type=float):
            return type(entryCols[catDict[fieldName]])
        entryID = catalog('id',int)

        # position on the sky in degrees
        RA = catalog('ra')
        DEC = catalog('dec')
        
        # skip sources outside our margins
        if RA < RAmin-margin or RA > RAmax+margin or DEC < DECmin-margin or DEC > DECmax+margin:
            continue
        
        # Calculate the offsets of this source from our image's bottom left corner in pixels
        # (which might be negative or byeond our image bounds because of the margins)
        xoffset = (RA - RAmin)*deg2arcsec/args.pixel_scale*RAscale
        yoffset = (DEC - DECmin)*deg2arcsec/args.pixel_scale

        # Look up redshift
        z = catalog('redshift')
        
        # Look up source AB magnitude in the requested band
        abMag = catalog(args.band + '_ab')
        # Correct for extinction
        abMag += args.extinction*(args.airmass - 1)
        # Calculate total detected signal in electrons
        flux = args.exposure_time*args.zero_point*math.pow(10,-0.4*(abMag-24))
        # Skip objects whose total flux is below our pixel threshold
        if flux < pixelCut:
            continue
        
        # Look up the component flux relative normalizations
        diskFluxNorm = catalog('fluxnorm_disk')
        bulgeFluxNorm = catalog('fluxnorm_bulge')
        agnFluxNorm = catalog('fluxnorm_agn')
        totalFluxNorm = diskFluxNorm + bulgeFluxNorm + agnFluxNorm

        # Calculate the disk and bulge fluxes to simulate
        if args.no_disk:
            diskFlux = 0
        else:
            diskFlux = flux*diskFluxNorm/totalFluxNorm
        if args.no_bulge:
            bulgeFlux = 0
        else:
            bulgeFlux = flux*bulgeFluxNorm/totalFluxNorm
        if diskFlux == 0 and bulgeFlux == 0:
            continue
        
        # Get disk component parameters
        if diskFlux > 0:
            hlr_d = catalog('DiskHalfLightRadius') # in arcsecs
            pa_d = catalog('pa_disk') # position angle in degrees
            a_d = catalog('a_d') # major axis length in arcsecs
            b_d = catalog('b_d') # minor axis length in arcsecs
            # Calculate sheared ellipse aspect ratio
            q_d = b_d/a_d # between 0.2 and 1
            # Convert position angle from degrees to radians
            pa_d = pa_d*deg2rad
            # Calculate bounding box in arcsecs without psf or pixel convolution
            (w_d,h_d) = sersicBounds(1,diskFlux+bulgeFlux,hlr_d,q_d,pa_d,sbCut)
        else:
            (w_d,h_d) = (0,0)
        
        # Get bulge component parameters
        if bulgeFlux > 0:
            hlr_b = catalog('BulgeHalfLightRadius') # in arcsecs
            pa_b = catalog('pa_bulge') # position angle in degrees
            a_b = catalog('a_b') # major axis length in arcsecs
            b_b = catalog('b_b') # minor axis length in arcsecs
            # Calculate sheared ellipse aspect ratio
            q_b = b_b/a_b # between 0.2 and 1
            # Convert position angle from degrees to radians
            pa_b = pa_b*deg2rad
            # Calculate bounding box in arcsecs without psf or pixel convolution
            (w_b,h_b) = sersicBounds(4,diskFlux+bulgeFlux,hlr_b,q_b,pa_b,sbCut)
        else:
            (w_b,h_b) = (0,0)

        # If a component is missing, set its nominal size and shape from the other component.
        if diskFlux == 0:
            (hlr_d,q_d,pa_d) = (hlr_b,q_b,pa_b)
        if bulgeFlux == 0:
            (hlr_b,q_b,pa_b) = (hlr_d,q_d,pa_d)

        # Combine the bulge and disk ellipticities
        (size,e1,e2) = combineEllipticities(hlr_d,q_d,pa_d,hlr_b,q_b,pa_b,bulgeFlux/(bulgeFlux+diskFlux))

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
        logger.info('Rendering input catalog line %d (entry id %d) with w x h = %d x %d' %
            (lineno,entryID,2*xhalf+1,2*yhalf+1))

        # Calculate the pixel coordinates of the stamp center.
        xstamp = 0.5*(bbox.xmin + bbox.xmax)
        ystamp = 0.5*(bbox.ymin + bbox.ymax)

        # Calculate the subpixel shift in arcsecs (not pixels!) of the source center
        # relative to the stamp center. Note that the resulting shift may be more than
        # one pixel in either direction because of the clipping operation above.
        xshift = (xoffset - (xstamp-0.5))*args.pixel_scale
        yshift = (yoffset - (ystamp-0.5))*args.pixel_scale

        if args.verbose:
            logger.info('    flux: %.3g electrons (%s-band AB %.1f)' % (flux,args.band,abMag))
            logger.info('  bounds: [%d:%d,%d:%d] pixels' % (bbox.xmin,bbox.xmax,bbox.ymin,bbox.ymax))
            logger.info('   shift: (%f,%f) arcsecs = (%f,%f) pixels' %
                (xshift,yshift,xshift/args.pixel_scale,yshift/args.pixel_scale))
            logger.info('    disk: frac = %f, hlr = %f arcsec, q = %f, beta = %f rad' %
                (diskFlux/flux,hlr_d,q_d,pa_d))
            logger.info('   bulge: frac = %f, hlr = %f arcsec, q = %f, beta = %f rad' %
                (bulgeFlux/flux,hlr_b,q_b,pa_b))
            logger.info('     agn: frac = %f' % (agnFluxNorm/flux))
            logger.info('    bbox: disk (%.1f,%.1f) bulge (%.1f,%.1f) psf %.1f arcsec' %
                (w_d,h_d,w_b,h_b,psfSize))
            logger.info('    size: %.2f pixels' % size)
            logger.info('   shear: (g1,g2) = (%.6f,%.6f)' % (e1,e2))
        
        # Define the nominal source parameters for rendering this object within its stamp
        params = {
            'total_flux': diskFlux + bulgeFlux,
            'f_d': diskFlux/(diskFlux+bulgeFlux), 'f_b': bulgeFlux/(diskFlux+bulgeFlux), 
            'x_d': xshift, 'y_d': yshift, 'hlr_d': hlr_d, 'q_d': q_d, 'beta_d': pa_d,
            'x_b': xshift, 'y_b': yshift, 'hlr_b': hlr_b, 'q_b': q_b, 'beta_b': pa_b,
            'dx': 0., 'dy': 0., 'relsize': 1., 'dbeta': 0.,
            'g1': args.g1, 'g2': args.g2
        }

        # Render the nominal stamps for this galaxy
        gal = createSource(**params)
        nominal = createStamp(gal,psf,pix,bbox)

        # Create a mask for pixels above threshold
        mask = createMask(nominal,pixelCut,args)
        if mask.array.sum() == 0:
            # this stamp has no pixels above threshold
            logger.info('*** line %d (id %d) is below threshold' % (lineno,entryID))
            continue
        trimmed = mask.bounds
        if args.verbose:
            logger.info(' trimmed: [%d:%d,%d:%d] pixels' %
                (trimmed.xmin,trimmed.xmax,trimmed.ymin,trimmed.ymax))

        # Add the nominal galaxy to the full field image after applying the threshold mask
        # (the mask must be the second term in the product so that the result is double precision)
        maskedNominal = nominal[trimmed]*mask
        fieldOverlap = trimmed & field.bounds
        if fieldOverlap.area() == 0:
             # this stamp's mask falls completely outside our field
            logger.info('*** line %d (id %d) does not overlap field' % (lineno,entryID))
            continue
        field[fieldOverlap] += maskedNominal[fieldOverlap]

        # Remember the nominal stamp (clipped to the field) for overlap calculations.
        stampList.append(maskedNominal[fieldOverlap])
        
        # Calculate this object's nominal flux S/N ratio at full depth using only masked pixels.
        # Note that this value cannot be reproduced from the saved stamp when a stamp is clipped
        # to the field boundary (use --no-clip to disable this).
        snr = signalToNoiseRatio(maskedNominal,args.exposure_time*skyRate)
        if args.verbose:
            logger.info('     S/N: %.6f' % snr)

        # Initialize the datacube of stamps that we will save for this object
        datacube = [ ]
        partialsArray = [ ]
        # Save the nominal (masked and trimmed) stamp
        assert saveStamp(datacube,maskedNominal,args)

        if args.partials:
            # Calculate the denominator array for our Fisher matrix elements
            fisherDenominator = maskedNominal[fieldOverlap].array + args.exposure_time*skyRate
            # Specify the amount to vary each parameter for partial derivatives
            # (we don't use a dictionary here since we want to control the order)
            variations = [
                ('f_d',0.01), ('f_b',0.01),
                ('dx',args.pixel_scale/3.),('dy',args.pixel_scale/3.),
                ('relsize',0.05),
                ('g1',0.03), ('g2',0.03)
            ]
            # the shape measurement parameters must always be the last 2 variations
            # since we make this assumption when slicing the covariance matrix below
            assert variations[-2][0] == 'g1'
            assert variations[-1][0] == 'g2'
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
                # append this partial to our datacube after trimming and masking
                maskedPartial = partial[trimmed]*mask
                assert saveStamp(datacube,maskedPartial,args)
                # remember this partial's numpy image array
                partialsArray.append(maskedPartial[fieldOverlap].array)
            # calculate the Fisher matrix images for this object
            nvar = len(partialsArray)
            nfisher = ((nvar+1)*nvar)/2
            (h,w) = partialsArray[0].shape
            fisherImage = numpy.zeros((nfisher,h,w))
            index = 0
            for i in range(nvar):
                for j in range(i,nvar):
                    fisherImage[index] = partialsArray[i]*partialsArray[j]/fisherDenominator
                    index += 1
            fisherImagesList.append(fisherImage)

        # Add a new HDU with a datacube for this object's stamps
        # We don't use compression = 'gzip_tile' for now since it is lossy
        # and mathematica cannot Import it.
        galsim.fits.writeCube(datacube, hdu_list = hduList)

        # Add a catalog entry for this galaxy
        entry = [entryID,xoffset,yoffset,abMag,flux/args.exposure_time,size,e1,e2,
            bulgeFlux/(diskFlux+bulgeFlux),z,snr]
        outputCatalog.append(entry)

        nkeep += 1
        logger.info("saved entry id %d as stamp %d" % (entryID,nkeep))

    # Loop over all saved objects to test for overlaps and build overlap groups
    (groupID,groupSize) = analyzeOverlaps(stampList)

    # Add group id to output catalog
    for (i,entry) in enumerate(outputCatalog):
        entry.append(groupID[i])

    # Do shape measurement error analysis for each galaxy
    purities = (0,0.5,0.9)
    regionsList = [ ]
    for i in range(nkeep):
        (errors,regions) = shapeErrorsAnalysis(nvar,stampList[i],fisherImagesList[i],field,purities)
        outputCatalog[i].extend(errors)
        regionsList.append(regions)

    # Save the regions for each object
    outname = args.output + '_regions.fits'
    logger.info('Saving regions to %r' % outname)
    galsim.fits.writeMulti(regionsList,outname)

    # Save group sizes to a file
    outname = args.output + '_groups.dat'
    out = open(outname,'w')
    ntot = 0
    for (i,n) in enumerate(groupSize):
        if n > 0:
            print >>out,i,n
            ntot += n
    out.close()
    assert ntot == len(groupSize)

    # Write the full field image without noise
    if args.save_field:
        # First without noise
        outname = args.output + '_field.fits'
        logger.info('Saving full field to %r' % outname)
        galsim.fits.write(field,outname)

    # Write the full field image with random noise added
    if args.save_noise:
        rng = galsim.BaseDeviate(123)
        noise = galsim.PoissonNoise(rng,sky_level = args.exposure_time*skyRate)
        field.addNoise(noise)
        outname = args.output + '_noise.fits'
        logger.info('Saving full field with noise added to %r' % outname)
        galsim.fits.write(field,outname)

    # Write the object stamp datacubes
    if args.stamps:
        outname = args.output + '_stamps.fits'
        logger.info('Saving %d stamps to %r' % (nkeep,outname))
        galsim.fits.writeFile(outname, hduList)

    # Close the input catalog
    cat.close()

    # Write the output catalog from memory
    outname = args.output + '_catalog.dat'
    out = open(outname,'w')
    for entry in outputCatalog:
        print >>out, ' '.join(map(str,entry))
    out.close()
    logger.info('Wrote %d of %d catalog entries to %r' % (nkeep,lineno,outname))

if __name__ == "__main__":
    main()
