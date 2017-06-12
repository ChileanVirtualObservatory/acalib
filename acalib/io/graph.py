from astropy import log
import numpy as np
from astropy.wcs import wcs
#from mayavi import mlab
from astropy.nddata import support_nddata
from ..core.analysis import rms
import matplotlib.pyplot as plt

#TODO: complete the nddata support (i.e. data, meta...)

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
    elif data.ndim == 3:
        if contour:
            visualize_contour3D(data,wcs,unit)
        else:
            visualize_volume(data,wcs,unit)
    else:
        log.error("Data dimensions must be between 1 and 3")

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

# TODO: Remove hardocded stuff
@support_nddata
def visualize_volume(data,wcs=None,unit=None):
     if wcs is None:
        log.error("WCS is needed by this function")
     figure = mlab.figure('Volume Plot')
     mesh=get_mesh(data)
     xi,yi,zi=mesh
     ranges=axes_ranges(data,wcs)
     grid = mlab.pipeline.scalar_field(xi, yi, zi, data)
     mmin = data.min()
     mmax = data.max()
     mlab.pipeline.volume(grid)#,vmin=mmin, vmax=mmin)
     ax=mlab.axes(xlabel="VEL [km/s] ",ylabel="DEC [deg]",zlabel="RA [deg]",ranges=ranges,nb_labels=5)
     ax.axes.label_format='%.3f'
     mlab.colorbar(title='flux', orientation='vertical', nb_labels=5)
     mlab.show()

@support_nddata
def visualize_contour3D(data,wcs=None,unit=None):
     if wcs is None:
        log.error("WCS is needed by this function")
     figure = mlab.figure('Contour Plot')
     mesh=get_mesh(data)
     xi,yi,zi=mesh
     ranges=axes_ranges(data,wcs)
     mmin = data.min()
     mmax = data.max()
     mlab.contour3d(xi,yi,zi,data,transparent=True,contours=10,opacity=0.5)
     ax=mlab.axes(xlabel="VEL [km/s] ",ylabel="DEC [deg]",zlabel="RA [deg]",ranges=ranges,nb_labels=5)
     ax.axes.label_format='%.3f'
     mlab.colorbar(title='flux', orientation='vertical', nb_labels=5)
     mlab.show()


def plot_snr_estimation(target,snr_results):
    fig = plt.figure(figsize=(10,4))
    # Unpack results
    (slimit,rms,x,y,v,n,p)=snr_results
    # Plot
    ax = fig.add_subplot(1,1,1)
    ax.plot(x,y,color='b')
    ax.plot(x,x,color='b',linestyle="--")
    ax.axvline(x=slimit,color='r',label="$"+str(slimit)+" \sigma$")
    ax.legend(loc=7)
    ax.set_title(target)
    ax.set_xlabel("SNR ($ \\tau / \sigma$)")
    ax.set_ylabel('RMS', color='b')
    for tl in ax.get_yticklabels():
        tl.set_color('b')
    axp = ax.twinx()
    axp.plot(x[:-1],v,color='grey')
    axp.set_ylabel('$\Delta$RMS', color='grey')
    for tl in axp.get_yticklabels():
        tl.set_color('grey')

