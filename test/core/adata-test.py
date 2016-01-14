import sys
sys.path.append('../../')

from acalib.core import *
from os import walk
from astropy import log
from matplotlib import pyplot as plt


binpath='../../bindata/fits/cubes/'
fn=binpath+'M100line.image.fits'

# Load from container
c = AContainer()
log.info("Loading "+fn)
c.load_from_fits(fn)

dt=c.adata[0]


# Simple Checks


shape=dt.shape()
log.info('Shape = '+str(shape))

pixels=dt.count()
log.info('Pixels = '+str(pixels))

tflux=dt.flux()
log.info('Total Flux =  '+str(tflux))

vflux=dt.variance()
log.info('Flux Variance =  '+str(vflux))

mmax=dt.max()
mmin=dt.min()
log.info('Max/Min = '+str(mmax)+', '+str(mmin))

rms=dt.estimate_rms()
log.info('RMS =  '+str(rms))

# TODO: Scale and rotate should be fixed and tested.

log.info('WCS axes = '+str(dt.axis_names()))

#WCS Limits
l1=dt.wcs_limits(0)
l2=dt.wcs_limits(1)
l3=dt.wcs_limits(2)
log.info('WCS limits = '+str(l1)+' '+str(l2)+' '+str(l3)+' ')

iw=dt.index_to_wcs((20,300,300))
log.info('WCS center= '+str(iw))

# We need to try these
#dt.index_from_window(self,wcs_center,wcs_window):
#dt.index_features(self,lower=None,upper=None):
#dt.features(self,lower=None,upper=None):

plt.figure(1)
plt.subplot(1,3,1)
plt.imshow(dt.stack())
plt.subplot(1,3,2)
plt.imshow(dt.stack(axis=1))
plt.subplot(1,3,3)
plt.imshow(dt.stack(axis=2))

(y_min,y_fact)=dt.standarize()
tflux=dt.flux()
log.info('Normalized Total Flux =  '+str(tflux))
dt.unstandarize((y_min,y_fact)

plt.show()

log.info(str(dt.cut(lower=(20,300,300),upper=(22,302,302))))
