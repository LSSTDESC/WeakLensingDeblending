import numpy as np 
import os
endpoint1 = (-.5 + -.25)/2
endpoint2 = -endpoint1
job_number = 2
for i,x in enumerate(np.linspace(endpoint1,endpoint2, 4)):
    for j,y in enumerate(np.linspace(endpoint1,endpoint2, 4)):
        cmd = './simulate.py --catalog-name OneDegSq.fits --survey-name LSST --image-width 4500 --image-height 4500 --output-name section{0}{1}{2} --ra-center {3} --dec-center {4} --calculate_bias --cosmic-shear-g1 0.01 --cosmic-shear-g2 0.01 --verbose'.format(i,j,job_number,x,y)
        slac_cmd = 'bsub -W 10:00 -o "output{0}{1}{2}.txt" -r "{3}"'.format(i,j,job_number,cmd)
        os.system(slac_cmd)