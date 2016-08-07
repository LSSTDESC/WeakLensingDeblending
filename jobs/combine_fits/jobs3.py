import numpy as np 
import subprocess
import os
import matplotlib.pyplot as plt
import galsim
import copy 
from astropy.table import Table
import fitsio

#directories that would be using 
WLD = '/Users/Ismael/code/lensing/WeakLensingDeblending/'
repo = '/Users/Ismael/code/lensing/repo/'
AEGIS = '/Users/Ismael/aegis/WeakLensingDeblending/'
intermediate_fits = '/Users/Ismael/aegis/WeakLensingDeblending/data/intermediate_fits/'
intermediate_fits_slac = '/nfs/slac/g/ki/ki19/deuce/AEGIS/ismael/WeakLensingDeblending/data/intermediate_fits'
os.chdir(intermediate_fits_slac)

#also have to write a new header to the fits file. 
job_number = 3
pixel_scale = .2

#read each. for given job_number. this works for the 1/16 setup. 
#ra: -.5 - .5 left to right, 
#dec: -.5 - .5 down to up, 
endpoint1 = (-.5 + -.25)/2
endpoint2 = -endpoint1
tables = []
#also have to combine the full galaxy postage stamp. 
stamps,past_i,past_j = None,None,None
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

#save full stamp 
from astropy.io import fits
new_fits_name_image = 'SimOSD{0}_image.fits'.format(job_number)
f = fits.PrimaryHDU(full_stamp)
f.writeto(new_fits_name_image)

from astropy.table import vstack 
Table = vstack(tables)

new_fits_name_table = 'SimOSD{0}_table.fits'.format(job_number)
Table.write(new_fits_name_table)

import astropy.io.fits as fits
#have to adjust to a correct header. 
f = fits.open(new_fits_name_table)
f_sample = fits.open('section00{0}.fits'.format(job_number))  #sample section of the job_number. 
f[0].header = f_sample[0].header
f[0].header['E_HEIGHT'] = 18000
f[0].header['GE_WIDTH'] = 18000


subprocess.call('rm {0}'.format(new_fits_name_table), shell=True) #delete older one so no problems at overwriting. 
f.writeto(new_fits_name_table)