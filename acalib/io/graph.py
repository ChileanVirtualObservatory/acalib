from astropy import log
import numpy as np
from astropy.wcs import wcs
from astropy.nddata import support_nddata
from ..core.analysis import rms
import matplotlib.pyplot as plt

#TODO: complete the nddata support (i.e. data, meta...)
#TODO: make animation possible again

@support_nddata
def visualize(data,wcs=None,unit=None,contour=False):
    """
    Generic function to visualize data, line-plot for 1D and image for 2D.

    Parameters
    ----------
    data : numpy.ndarray or astropy.nddata.NDData or astropy.nddata.NDDataRef
        Astronomical image

    wcs : astropy.wcs.WCS
        World Coordinate System from the image (not needed if contained in NDData)

    unit : astropy.unit
        Image units (not needed if contained in NDData)

    contour : numpy.ndarray
        For plotting Contourns
    """
    if data.ndim == 1:
        return visualize_plot(data,wcs,unit)
    elif data.ndim == 2:
        return visualize_image(data,wcs,unit,contour)
    else:
        log.error("Data dimensions must be between 1 and 2")

@support_nddata
def visualize_plot(data,wcs=None,unit=None):
    """
    Plot 1D data for astronomical data.

    Parameters
    ----------
    data : numpy.ndarray or astropy.nddata.NDData
        Astronomical image

    wcs : astropy.wcs.WCS
        World Coordinate System from the image (not needed if contained in NDData)

    unit : astropy.unit
        Image units (not needed if contained in NDData)
    """
    if wcs is None:
        plt.plot(data)
        plt.ylabel(unit)
    else:
        #TODO: Implement x vector, but check why the wcs cannot be onedimensional!
        plt.plot(data)
        plt.ylabel(unit)
        plt.xlabel(wcs.axis_type_names[0])
    #plt.show()

@support_nddata
def visualize_image(data,wcs=None,unit=None,contour=False):
    """
    Plot 2D astronimical data.

    Parameters
    ----------
    data : numpy.ndarray or astropy.nddata.NDData or astropy.nddata.NDDataRef
        Astronomical image

    wcs : astropy.wcs.WCS
        World Coordinate System from the image (not needed if contained in NDData)

    unit : astropy.unit
        Image units (not needed if contained in NDData)

    contour : numpy.ndarray
        For plotting Contourns
    """
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
        arms=rms(data)
        dmax=data.max()
        crs=np.arange(1,dmax/arms)
        plt.contour(data,levels=arms*crs,alpha=0.5)
    plt.show()
