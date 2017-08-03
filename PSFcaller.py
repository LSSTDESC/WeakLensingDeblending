def PSF_at_position():
    ''' I'm drafting this function to call the PSF at a certain position, and input it as a
    Survey.psf_model object into WeakLensingDeblending. The function should return a GSObject '''

     from angles import r2d
     import sys, os, matplotlib,  galsim
     import numpy as np
     from lsst.daf.persistence import Butler
     from lsst.afw.geom import Point2D
     matplotlib.use("pdf")
     import matplotlib.pyplot as plt
     import seaborn as sns;sns.set_style('darkgrid')

     if len(sys.argv)<2:
         sys.stderr.write("Syntax: python test-psf.py  repo_path\n")
         exit()

     repo_rel_path = sys.argv[1]
     repo_abs_path = os.path.abspath(repo_rel_path)

     if not os.path.exists(repo_abs_path):
         sys.stderr.write("Nothing found at {}\n".format(repo_abs_path))
         exit()

     # Create a data butler which provides access to a data repository.
     butler = Butler(repo_abs_path)
     print "Butler summoned (i.e. we have loaded the data repository)", repo_abs_path
     ccd_exposures = butler.queryMetadata("calexp", format=['visit', 'ccd'])
     ccd_exposures = [(visit,ccd) for (visit,ccd) in ccd_exposures if butler.datasetExists("calexp", visit=visit, ccd=ccd)]

     print "Found {} (exposure,ccd) pairs".format(len(ccd_exposures))
     for ccd in range(104):
         if ccd ==1:
             old_X, old_Y = X, Y
         visit, ccd = ccd_exposures[ccd]
         calexp = butler.get("calexp", visit=visit, ccd=ccd, immediate=True)
         detector = calexp.getDetector()
         print 'visit', visit, 'ccd', ccd, detector.getName()
         width, height, psf, wcs = calexp.getWidth(), calexp.getHeight(), calexp.getPsf(), calexp.getWcs()
         nx = 20
         ny = int(nx * (height/width))
         X,Y= (np.zeros((ny,nx)) for _ in range(2))
         for i,x in enumerate(np.linspace(0, width, nx)):
             for j,y in enumerate(np.linspace(0, height, ny)):
                 point = Point2D(x,y)
                 image = psf.computeKernelImage(point)
                 galsim_image = galsim.Image(image.getArray())
                 shape_data = galsim.hsm.FindAdaptiveMom(galsim_image,strict=False)
                 if shape_data.error_message=='':
                     pass
                 else:
                     print shape_data.error_message
                     continue
                 FpPos = wcs.pixelToSky(point)
                 X[j,i], Y[j,i] = FpPos.getLongitude(), FpPos.getLatitude()
         if ccd == 0:
             pass
         else:
             old_X, old_Y = np.append(X,old_X), np.append(Y,old_Y)
