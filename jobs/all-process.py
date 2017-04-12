import numpy as np 
import os
import sys 
import argparse


def main():

    parser = argparse.ArgumentParser(description=('Simulate 16 different regions from a square' 
                                                  'degree and analyze their combination in'
                                                  'SExtractor.'),
                                     formatter_class=(
                                     argparse.ArgumentDefaultsHelpFormatter))

    parser.add_argument('--simulate', action='store_true',
                        help=('Simulates the requested 16 regions with given job_number'))

    parser.add_argument('--combine', action='store_true',
                        help=('Combines 16 regions with given job_number'))

    parser.add_argument('--extract', action='store_true',
                        help=('Classifies into detected and ambiously blended.'))

    parser.add_argument('--job-number', required=True,
                        type=int,
                        help=('Job id number'))

    parser.add_argument('--cosmic-shear-g1', required=True,
                        type=float,
                        help=('cosmic shear g1'))

    parser.add_argument('--cosmic-shear-g2', required=True,
                        type=float,
                        help=('cosmic shear g2'))

    parser.add_argument('--noise-seed', required=True,
                        type=int,
                        help=('noise seed to use when adding noise to image SExtracted.'))

    args = parser.parse_args()

    ###################################################################
    #some constants used throughout
    pixel_scale = .2
    total_width = 18000
    total_height = 18000


    endpoint1 = (-.5 + -.25)/2
    endpoint2 = -endpoint1

    ###################################################################
    #names to be used. 

    inputs = dict(
    table = '/Users/Ismael/aegis/data/intermediate_fits/SimOSD{0}_table.fits'.format(args.job_number),
    noise_image = '/Users/Ismael/aegis/data/intermediate_fits/SimOSD{0}_image.fits'.format(args.job_number),
    output_detected = '/Users/Ismael/aegis/data/sextractor_runs/out{0}.cat'.format(args.job_number),
    final_fits = '/Users/Ismael/aegis/data/fits_files/final_fits{0}.fits'.format(args.job_number),
    config_file = '/Users/Ismael/aegis/data/sextractor_runs/default.sex',
    WLD = '/Users/Ismael/aegis/WeakLensingDeblending/',
    sample_fits = '/Users/Ismael/aegis/data/intermediate_fits/section00{0}.fits'.format(args.job_number),
    intermediate_fits = '/Users/Ismael/aegis/data/intermediate_fits', 
    )

    inputs_slac = dict()
    for f in inputs: 
        l = inputs[f].split("/")
        index_aegis = l.index("aegis")
        str_slac = "/".join(l[index_aegis+1:])
        slac_file = '{0}{1}'.format('/nfs/slac/g/ki/ki19/deuce/AEGIS/ismael/',str_slac)
        inputs_slac[f] = slac_file


    ######################################################################################################################################
    #simulate the 16 regions. 
    os.chdir(inputs_slac['WLD'])

    if args.simulate:
        for i,x in enumerate(np.linspace(endpoint1,endpoint2, 4)):
            for j,y in enumerate(np.linspace(endpoint1,endpoint2, 4)):
                cmd = './simulate.py --catalog-name OneDegSq.fits --survey-name LSST --image-width 4500 --image-height 4500 --output-name section{0}{1}{2} --ra-center {3} --dec-center {4} --calculate_bias --cosmic-shear-g1 {5} --cosmic-shear-g2 {6} --verbose'.format(i,j,args.job_number,x,y,args.cosmic_shear_g1,args.cosmic_shear_g2)
                slac_cmd = 'bsub -W 10:00 -o "output{0}{1}{2}.txt" -r "{3}"'.format(i,j,args.job_number,cmd)
                os.system(slac_cmd)



    ######################################################################################################################################
    #combine 16 regions into image and table - add noise to the image. 

    if args.combine

        os.chdir(inputs_slac['intermediate_fits'])

        #read each. for given job_number. this works ONLY for the 1/16 setup. 
        #ra: -.5 - .5 left to right, 
        #dec: -.5 - .5 down to up, 

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
                table['dx']+=x*total_width*pixel_scale #use lsst pixel scale 
                table['dy']+=y*total_height*pixel_scale
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
        reader = descwl.output.Reader(inputs_slac['sample_fits'])
        results = reader.results
        generator = galsim.random.BaseDeviate(seed = noise_seed)
        noise = galsim.PoissonNoise(rng = generator, sky_level = results.survey.mean_sky_level)
        full_stamp_galsim = galsim.Image(array=full_stamp,wcs=galsim.PixelScale(pixel_scale),bounds=galsim.BoundsI(xmin=0, xmax=total_width-1, ymin=0, ymax=total_height-1)
        full_stamp_galsim.addNoise(noise)


        #save full stamp 
        import astropy.io.fits as fits
        f = fits.PrimaryHDU(full_stamp_galsim.array)
        f.writeto(inputs_slac['noise_image'])


        from astropy.table import vstack 
        Table = vstack(tables)
        Table.write(inputs_slac['table'])

        #have to adjust to a correct header. 
        f = fits.open(inputs_slac['table'])
        f_sample = fits.open(inputs_slac['sample_fits'])  #sample section of the job_number. 
        f[0].header = f_sample[0].header
        f[0].header['E_HEIGHT'] = total_height
        f[0].header['GE_WIDTH'] = total_width
        # f[0].header['NAXIS1'] = 18000
        # f[0].header['NAXIS2'] = 18000
        subprocess.call('rm {0}'.format(inputs_slac['table']), shell=True) #delete older one so no problems at overwriting. 
        f.writeto(inputs_slac['table'])


    ######################################################################################################################################
    #source extract and detect ambiguous blends 

    if args.extract:

        cmd = 'sex {0} -c {1} -CATALOG_NAME {2}'.format(inputs_slac['noise_image'],inputs_slac['config_file'],inputs_slac['output_detected'])
        os.system(cmd)

        cat = descwl.output.Reader(inputs_slac['table']).results
        table = cat.table
        detected,matched,indices,distance = cat.match_sextractor(inputs_slac['output_detected'])

        #convert to arcsecs and relative to image_center 
        detected['X_IMAGE'] = (detected['X_IMAGE'] - 0.5*total_width - 0.5)*pixel_scale
        detected['Y_IMAGE'] = (detected['Y_IMAGE'] - 0.5*total_height - 0.5)*pixel_scale

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
        subprocess.call('cp {0} {1}'.format(inputs_slac['table'],inputs_slac['final_fits']),shell=True)
        f = fits.open(inputs_slac['final_fits'])

        #replace table and write file
        f[1]= astropy.io.fits.table_to_hdu(table)
        subprocess.call('rm {0}'.format(inputs_slac['final_fits']),shell=True)
        f.writeto(inputs_slac['final_fits']) #overwrites the existing one. 


if __name__=='__main__':
    main()