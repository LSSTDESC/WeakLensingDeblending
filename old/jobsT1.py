import numpy as np 
import os

job_number = 'T1'
for i in range(5):
    cmd = './simulate.py --catalog-name OneDegSq.fits --survey-name LSST --image-width 18000 --image-height 18000 --output-name section_final{0}{1} --calculate_bias --cosmic-shear-g1 0.01 --cosmic-shear-g2 0.01 --verbose'.format(i,job_number)
    slac_cmd = 'bsub -W 100:00 -o "output_all{0}{1}.txt" -r "{2}"'.format(i,job_number,cmd)
    os.system(slac_cmd)