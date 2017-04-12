import numpy as np 
import subprocess
import os
import matplotlib.pyplot as plt
import galsim
import copy 
from astropy.table import Table
import fitsio
import descwl

#directories that would be using 
WLD = '/Users/Ismael/code/lensing/WeakLensingDeblending/'
repo = '/Users/Ismael/code/lensing/repo/'
AEGIS = '/Users/Ismael/aegis/'
intermediate_fits = '/Users/Ismael/aegis/data/intermediate_fits/'
aegis_slac = '/nfs/slac/g/ki/ki19/deuce/AEGIS/ismael/'
intermediate_fits_slac = '/nfs/slac/g/ki/ki19/deuce/AEGIS/ismael/data/intermediate_fits'
os.chdir(intermediate_fits_slac)

#also have to write a new header to the fits file. 
job_number = 2
pixel_scale = .2
sample_fits_name = 'section00{0}.fits'.format(job_number) #use to get some info. 
noise_seed = 0 #when adding noise to the image. 


#read each. for given job_number. this works for the 1/16 setup. 
#ra: -.5 - .5 left to right, 
#dec: -.5 - .5 down to up, 
endpoint1 = (-.5 + -.25)/2
endpoint2 = -endpoint1
tables = []
#also have to combine the full galaxy postage stamp. 
stamps,past_i = None,None
for i,x in enumerate(np.linspace(endpoint1,endpoint2, 4)):
    for j,y in enumerate(np.linspace(endpoint1,endpoint2, 4)):
        file_name = 'section{0}{1}{2}.fits'.format(i,j,job_number)
        table = Table.read(file_name)
        fits_section = fitsio.FITS(file_name)
        stamp = fits_section[0].read()
        #adjust corresponding entries to the absolute image center 
        table['dx']+=x*18000*pixel_scale #use lsst pixel scale 
        table['dy']+=y*18000*pixel_scale
        tables.append(table)
        if stamps is None: 
            stamps = []
            stamps.append(stamp) 
            past_i = i 
        else: 
            if past_i == i: 
                stamps[i] = np.vstack([stamp,stamps[i]])
            else: 
                past_i = i 
                stamps.append(stamp)
#combine stamps list into final stamp 
full_stamp = np.hstack(stamps)


#take any noise_seed and add noise to the image generated
reader = descwl.output.Reader(sample_fits_name)
results = reader.results
generator = galsim.random.BaseDeviate(seed = noise_seed)
noise = galsim.PoissonNoise(rng = generator, sky_level = results.survey.mean_sky_level)
full_stamp_galsim = galsim.Image(array=full_stamp,wcs=galsim.PixelScale(pixel_scale),bounds=galsim.BoundsI(xmin=0, xmax=17999, ymin=0, ymax=17999))
full_stamp_galsim.addNoise(noise)


#save full stamp 
import astropy.io.fits as fits
new_fits_name_image = 'SimOSD{0}_image.fits'.format(job_number)
f = fits.PrimaryHDU(full_stamp_galsim.array)
f.writeto(new_fits_name_image)


# fits = fitsio.FITS('section00{0}.fits')
# fits[0].read() = full_stamp #replace postage stamp 

# #all galaxies, simulated survey image 
# print fits[0]
# img = fits[0][:,:] #notice to get 512x512 image we do this
# img = fits[0].read() #can also do this. 

# from astropy.io import fits


from astropy.table import vstack 
Table = vstack(tables)

new_fits_name_table = 'SimOSD{0}_table.fits'.format(job_number)
Table.write(new_fits_name_table)

#have to adjust to a correct header. 
f = fits.open(new_fits_name_table)
f_sample = fits.open(sample_fits_name)  #sample section of the job_number. 
f[0].header = f_sample[0].header
f[0].header['E_HEIGHT'] = 18000
f[0].header['GE_WIDTH'] = 18000
# f[0].header['NAXIS1'] = 18000
# f[0].header['NAXIS2'] = 18000
subprocess.call('rm {0}'.format(new_fits_name_table), shell=True) #delete older one so no problems at overwriting. 
f.writeto(new_fits_name_table)