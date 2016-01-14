import sys
sys.path.append('../../')

from acalib.core import *
from os import walk
from astropy import log

binpath='../../bindata/fits/cubes/'
fn=binpath+'M100line.image.fits'

# Load from container
c = AContainer()
log.info("Loading "+fn)
c.load_from_fits(fn)

dt=c.adata[0]


# Simple Data

rms=dt.estimate_rms()
log.info('RMS =  '+str(rms))

tflux=dt.flux()
log.info('Total Flux =  '+str(tflux))

pixels=dt.count()
log.info('Pixels = '+str(pixels))

shape=dt.shape()
log.info('Shape = '+str(shape))



#WCS Limits
l1=dt.get_wcs_limits(0)
l2=dt.get_wcs_limits(1)
l3=dt.get_wcs_limits(2)
log.info('WCS limits = '+str(l1)+' '+str(l2)+' '+str(l3)+' ')



#    def scale(self, scale):
#        dim = 0
#        start_time = time.time()
#        if (scale == 1):
#            return self.data
#        elif (scale < 1):
#            new_data = self.data[::1/scale, ::1/scale, ::1/scale]
#            return (new_data/np.sum(new_data))*np.sum(self.data)
#        else:
#            new_data = np.zeros((round(len(self.data)*scale),round(len(self.data[0])*scale), round(len(self.data[0][0])*scale)))
#            new_data[::scale, ::scale, ::scale] = self.data
#
#
#            new_data = interpolate(new_data, scale)
#
#            return new_data
#
#
#    def rotate(self, angle):
#        if (angle != 0):
#            new_data = self.data.rotate(angle, Image.BICUBIC,1)
#            return new_data
#        else:
#            return self.data


