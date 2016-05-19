import sys
sys.path.append('../../')

import acalib
import acalib.io.formats as io
import acalib.roi as ri
import matplotlib.pyplot as plt


#c = acalib.AContainer()
#io.load_to_cont('/home/cvalenzu/Escritorio/indexing_utfsm/ALMA00000085.fits',c)
#cube = c.primary
rd = ri.RoiDetect()


#print "Sketching Cube Spectra"
#spectra = rd.cube_spectra(cube,50)

#print "Plotting"
#frecs = range(950)
#plt.bar(frecs, spectra)
#plt.show()

#h1 = rd.vel_stacking(cube,0,30)

#h1t = h1.T
#plt.imshow(h1t, origin='lower')
#plt.show()

#rd.gaussian_mix(h1)
from astropy.io import fits


rd = ri.RoiDetect()
hdulist = fits.open('/home/cvalenzu/Escritorio/indexing_utfsm/drC-crop.fits')
ndd = hdulist[0]
data = ndd.data

fts = rd.gaussian_mix(data)

#import numpy as np
#from skimage.filter import threshold_adaptive
#image = data.astype(float)

#f = (image-np.min(image))/(np.max(image)-np.min(image))
#g = threshold_adaptive(f,360,method='mean',offset=0)

#plt.imshow(g, origin='lower')
#plt.show()