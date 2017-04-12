import numpy as np 
import subprocess
import os
import matplotlib.pyplot as plt
import galsim
import copy 
from astropy.table import Table
import fitsio
import argparse

#inputs, put all input files in a dict, 

#algorithm for ambiguous blends - returns rows of ambiguously blended objects in cat. 
#the difference is that it tags all objects <1. effective distance as blended. 


def main(): 

    parser = argparse.ArgumentParser(description=('Detects ambiguos blends from a given image and'
                                                  'catalog of image creates a new fits with table'
                                                  'containing which galaxies are ambiguosly'
                                                  'blended and detected. '),
                                 formatter_class=(
                                 argparse.ArgumentDefaultsHelpFormatter))

    parser.add_argument('--job-number', required=True,
                        type=int,
                        help=('Job id number'))

    args = parser.parse_args()

    #directories that would be using 
    WLD = '/Users/Ismael/code/lensing/WeakLensingDeblending/'
    WLFF = '/Users/Ismael/code/lensing/WLFF/'
    aegis = '/Users/Ismael/aegis/data/'
    SEx = '/Users/Ismael/aegis/data/sextractor_runs/'
    aegis_slac = '/nfs/slac/g/ki/ki19/deuce/AEGIS/ismael/data'
    temp_data = '/Users/Ismael/temp_data'
    os.chdir(WLD)
    import descwl

    inputs = dict(
    input_name_table = '/Users/Ismael/aegis/data/intermediate_fits/SimOSD{0}_table.fits'.format(args.job_number),
    input_name_image = '/Users/Ismael/aegis/data/intermediate_fits/SimOSD{0}_image.fits'.format(args.job_number),
    input_name_image_slac = '/nfs/slac/g/ki/ki19/deuce/AEGIS/ismael/data/intermediate_fits/SimOSD{0}_image.fits'.format(args.job_number),
    output_name = '/Users/Ismael/aegis/data/sextractor_runs/out{0}.cat'.format(args.job_number),
    config_file_slac = '/nfs/slac/g/ki/ki19/deuce/AEGIS/ismael/data/sextractor_runs/default.sex',
    final_fits = '/Users/Ismael/aegis/data/fits_files/final_fits{0}.fits'.format(args.job_number),
    output_name_slac='/nfs/slac/g/ki/ki19/deuce/AEGIS/ismael/data/sextractor_runs/out{0}.cat'.format(args.job_number),     
    )

    cmd = 'sex {0} -c {1} -CATALOG_NAME {2}'.format(inputs['input_name_image_slac'],inputs['config_file_slac'],inputs['output_name_slac'])
    print cmd

    input() 

    cat = descwl.output.Reader(inputs['input_name_table']).results
    table = cat.table
    detected,matched,indices,distance = cat.match_sextractor(inputs['output_name'])

    width, height = 18000,18000
    pixel_scale = .2 

    #convert to arcsecs and relative to image_center 
    detected['X_IMAGE'] = (detected['X_IMAGE'] - 0.5*width - 0.5)*pixel_scale
    detected['Y_IMAGE'] = (detected['Y_IMAGE'] - 0.5*height - 0.5)*pixel_scale

    #convert second moments arcsecs 
    detected['X2_IMAGE']*=pixel_scale**2 
    detected['Y2_IMAGE']*=pixel_scale**2 
    detected['XY_IMAGE']*=pixel_scale**2 


    # calculate size from moments X2_IMAGE,Y2_IMAGE,XY_IMAGE -> remember in pixel**2 so have to convert to arcsecs. 
    sigmas = []
    for x2,y2,xy in zip(detected['X2_IMAGE'],detected['Y2_IMAGE'],detected['XY_IMAGE']):
        second_moments = np.array([[x2,xy],[xy,y2]])
        sigma = np.linalg.det(second_moments)**(+1./4) #should be a plus 
        sigmas.append(sigma)
    SIGMA = Table.Column(name='SIGMA',data=sigmas)
    detected.add_column(SIGMA)
    # #add sizes of the corresponding entries from lsst to the detected table. 
    # SIGMA_M = Table.Column(name='SIGMA_M', data=table[indices]['sigma_m'])
    # detected.add_column(SIGMA_M) 

    ambg_blends = detected_ambiguous_blends(table, indices, detected)
    ambg_blends_indices = list(ambg_blends)
    ambg_blends_ids = list(table[ambg_blends_indices]['db_id'])


    #add columns to table of undetected and ambiguosly blended 
    ambigous_blend_column = []
    for i,gal_row in enumerate(table):
        if i in ambg_blends_indices:
            ambigous_blend_column.append(True)
        else: 
            ambigous_blend_column.append(False)
    column = Table.Column(name='ambig_blend',data=ambigous_blend_column)
    table.add_column(column)


    #create new copy of fits_file with this new table. 
    subprocess.call('cp {0} {1}'.format(inputs['input_name'],inputs['final_fits']),shell=True)
    f = fits.open(inputs['final_fits'])

    #replace table and write file
    f[1]= astropy.io.fits.table_to_hdu(table)
    subprocess.call('rm {0}'.format(inputs['final_fits']),shell=True)
    f.writeto(inputs['final_fits']) #overwrites the existing one. 

if __name__=='__main__':
    main()