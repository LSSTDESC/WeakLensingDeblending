import math

# pixel scale in arcsecs / pixel
PixelScale = 0.2

# size of square region to use in pixels, with its top-left corner at (ra,dec) = (0,0)
# (an LSST chip is 4096 pixels on a side)
FieldSize = 4096

# total signal per 15s exposure in ADU at AB magnitude 24.0
# use 711, 810 in LSST i,r-bands
fluxnorm = 711

# surface brightness isophote to cut at in ADU units per pixel where, in the case of an
# image stack, ADU is an average rather than a sum so that the sky decreases with stacking.
fcut = 2.0
f0 = fcut/(PixelScale*PixelScale)

# psf parameters to use when calculating bounding boxes
fwhm = 0.7 # arcsecs
beta = 3.0

# calculate r0 from fwhm
r0psf = 0.5*fwhm/math.sqrt(math.pow(2.,1./beta)-1)
cpsf = math.pi*f0*r0psf*r0psf/(beta-1.)

fin = open('gcat.dat','r')
fout = open('trim.dat','w')

twopi = 2*math.pi
deg2rad = math.pi/180.
deg2arcsec = 3600.

first = True
lineno = 0
count = 0
for line in open('gcat.dat'):
    lineno += 1
    cols = line.split()
    if first:
        names = cols
        first = False
    else:
        ra = float(cols[1])*deg2arcsec
        dec = FieldSize*PixelScale + float(cols[2])*deg2arcsec
        # ignore objects outside our window
        if ra > FieldSize*PixelScale or dec < 0:
            continue
        # process disk components
        hlr_d = float(cols[7]) # in arcsecs
        if hlr_d > 0:
            # Lookup extra params for this object
            r_ab = float(cols[22]) # AB magnitude in r band
            i_ab = float(cols[23]) # AB magnitude in i band
            pa_d = float(cols[9]) # position angle in degrees
            a_d = float(cols[17]) # major axis length in arcsecs
            b_d = float(cols[19]) # minor axis length in arcsecs
            # Calculate total flux in ADU units
            flux = fluxnorm*math.pow(10,24-i_ab)
            # Ignore objects whose total flux is below half of our per-pixel threshold
            if flux < 0.5*fcut:
                continue
            # Calculate sheared ellipse aspect ratio
            q_d = b_d/a_d # between 0.2 and 1
            # Convert position angle from degrees to radians
            pa_d = pa_d*deg2rad
            # Convert half-light radius of exponential profile to scale radius in arcsecs
            r_scale = hlr_d/1.67835
            # Calculate shear affine transform parameters
            g = (1-q_d)/(1+q_d)
            gp = g*math.cos(2*pa_d)
            gx = g*math.sin(2*pa_d)
            detM = 1 - gp*gp - gx*gx
            # Calculate the bounding box for the limiting isophote with no psf
            x = twopi*r_scale*r_scale*f0/(flux*detM)
            rcut = -r_scale*math.log(x)
            if rcut <= 0:
                # object is below our threshold without any psf
                continue
            dx = rcut*math.sqrt(((1+gp)*(1+gp)+gx*gx)/detM) # half width in arcsecs
            dy = rcut*math.sqrt(((1-gp)*(1-gp)+gx*gx)/detM) # half height in arcsecs
            # Calculate the psf padding
            arg = math.pow(cpsf/flux,-1./beta) - 1
            if arg > 0:
                rpad = r0psf*math.sqrt(arg)
            else:
                rpad = 0
            dx += rpad
            dy += rpad
            print >>fout,ra,dec,flux,hlr_d,q_d,pa_d,dx,dy,i_ab
            count += 1

print 'using',count,'of',lineno,'objects'

fout.close()
fin.close()
