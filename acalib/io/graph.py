from astropy import log
import numpy as np
from astropy.wcs import wcs
from astropy.nddata import support_nddata
from acalib import *
from acalib.core.indices import *
from acalib.core.utils import *


#TODO: complete the nddata support (i.e. data, meta...)
#TODO: make animation possible again

@support_nddata
def visualize(data,wcs=None,unit=None,contour=False):
    if data.ndim == 1:
        return visualize_plot(data,wcs,unit)
    elif data.ndim == 2:
        return visualize_image(data,wcs,unit,contour)
    else:
        log.error("Data dimensions must be between 1 and 2")

@support_nddata
def visualize_plot(data,wcs=None,unit=None):
    if wcs is None:
        plt.plot(data)
        plt.ylabel(unit)
    else:
        #TODO: Implement x vector, but check why the wcs cannot be onedimensional!
        plt.plot(data)
        plt.ylabel(unit)
        plt.xlabel(wcs.axis_type_names[0])
    plt.show()
         
@support_nddata
def visualize_image(data,wcs=None,unit=None,contour=False):
    if wcs is None:
        plt.imshow(data, origin='lower', cmap=plt.cm.gist_heat)
        cb=plt.colorbar()
        cb.ax.set_ylabel(unit)
    else:
        gax=plt.subplot(111,projection=wcs)
        plt.imshow(data, origin='lower', cmap=plt.cm.gist_heat)
        g0=gax.coords[0]
        g1=gax.coords[1]
        g0.set_axislabel(wcs.axis_type_names[0])
        g1.set_axislabel(wcs.axis_type_names[1])
        g0.grid(color='yellow', alpha=0.5, linestyle='solid')
        g1.grid(color='yellow', alpha=0.5, linestyle='solid')
        cb=plt.colorbar()
        cb.ax.set_ylabel(unit)
    if contour:
        rms=estimate_rms(data)
        dmax=data.max()
        crs=np.arange(1,dmax/rms)
        plt.contour(data,levels=rms*crs,alpha=0.5)
    plt.show()

