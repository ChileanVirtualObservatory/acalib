from astropy import log
import numpy as np
from astropy.wcs import wcs
#from mayavi import mlab
from astropy.nddata import support_nddata

from acalib import core, upi
from acalib.core import rms
from acalib.upi import axes, flux
import matplotlib.pyplot as plt
import ipyvolume.pylab as ipvlab


def _draw_spectra(data, wcs=None, unit=None,velocities=False):
    try:
        freq_axis = axes._get_axis(wcs, "FREQ")
    except ValueError:
        log.warning("Data does not have a spectral dimension")
        return
    ll = list(range(wcs.naxis))
    ll.remove(freq_axis)
    yvals = np.nansum(data, axis=tuple(ll))
    plt.ylabel(unit)
    if velocities == True:
        xvals = axes.spectral_velocities(data, wcs=wcs)
        plt.xlabel("VEL [Km/s]")
    else:
        dim = wcs.wcs.spec
        fqis = np.arange(data.shape[data.ndim - dim - 1])
        idx = np.zeros((fqis.size, data.ndim))
        idx[:, dim] = fqis
        vals = wcs.all_pix2world(idx, 0)
        xvals = vals[:, dim]
        plt.xlabel("FREQ [Hz]")
    plt.plot(xvals, yvals)


@support_nddata
def visualize_spectra(data, wcs=None, unit=None,velocities=False):
    if wcs is None:
        if data.ndim != 1:
            log.info("Only 1D data can be shown without WCS")
            return
        else:
            visualize_plot(data,wcs,unit)
            return
    _draw_spectra(data,wcs,unit,velocities)
    plt.show()

@support_nddata
def visualize(data,wcs=None,unit=None):
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
    """

    if data.ndim == 1:
        return visualize_plot(data,wcs,unit)
    elif data.ndim == 2:
        return visualize_image(data,wcs,unit)
    elif data.ndim == 3:
        if axes._is_spectra(data,wcs):
            visualize_spectra(data,wcs,unit)
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
    plt.show()



def _draw_image(data,gax,wcs,unit,cmap=plt.cm.gist_heat):
    plt.imshow(data, origin='lower', cmap=cmap)
    g0 = gax.coords[0]
    g1 = gax.coords[1]
    g0.set_axislabel(wcs.axis_type_names[0])
    g1.set_axislabel(wcs.axis_type_names[1])
    cb = plt.colorbar()
    cb.ax.set_ylabel(unit)


@support_nddata
def show_subcube(data,wcs=None,meta=None,mask=None,unit=None,lower=None,upper=None):
    plt.figure(figsize=(10,5))
    mslab=core.slab(data, lower, upper)
    #scube=reduction.cut(data,wcs=wcs,meta=meta,mask=mask,unit=unit,lower=lower,upper=upper)
    if data.ndim == 3:
        freq_axis = axes._get_axis(wcs, "FREQ")
        wcs=wcs.slice(mslab, numpy_order=True)
        wcsp=wcs.dropaxis(wcs.wcs.spec)
        im1=121
        #data2d= np.nansum(data, axis=(freq_axis))
        #scube2d=np.nansum(scube.data, axis=(freq_axis))
    else:
        wcs=wcs.slice(mslab, numpy_order=True)
        wcsp=wcs
        im1=111
        #im2=122
        #data2d=data
        #scube2d=scube.data
    gax=plt.subplot(im1,projection=wcsp)
    if data.ndim == 3:
        _draw_image(np.nansum(data[mslab],axis=freq_axis),gax,wcsp,unit,cmap=plt.cm.seismic)
    else:
        _draw_image(data[mslab], gax, wcsp, unit, cmap=plt.cm.seismic)
    #gax=plt.subplot(im2, projection=swcs2d)
    #_draw_image(scube2d, gax, swcs2d, unit, cmap=plt.cm.seismic)
    if (data.ndim == 3):
        #plt.subplot(223)
        #_draw_spectra(data,wcs=wcs,unit=unit)
        plt.subplot(122)
        _draw_spectra(data[mslab], wcs=wcs, unit=unit)
    plt.show()
    #return(scube)

def show_clumps(hrtree,node=0,print_numbers=True):
    fig = plt.figure(figsize=(10, 10))
    shape = hrtree.synth.data.shape
    parent=hrtree.find_parent(node,hrtree.tree)
    n_colors = hrtree.visible_clumps(parent,node)
    gs = gridspec.GridSpec(2, 2,width_ratios=[shape[2] / shape[0], 1], height_ratios=[shape[1] / shape[0], 1])
    tempXY = hrtree.synth.data.sum(axis=(0))
    nax = plt.subplot(gs[0])
    nax.imshow(tempXY, origin='lower', cmap=plt.cm.gray_r)
    nax.get_xaxis().set_visible(False)
    nax.get_yaxis().set_visible(False)

    tempXZ = hrtree.synth.data.sum(axis=(1))
    nax2 = plt.subplot(gs[2])
    nax2.imshow(tempXZ, origin='lower', cmap=plt.cm.gray_r)
    nax2.get_xaxis().set_visible(False)
    nax2.get_yaxis().set_visible(False)

    tempYZ = hrtree.synth.data.sum(axis=(2)).T
    nax3 = plt.subplot(gs[1])
    nax3.imshow(tempYZ, origin='lower', cmap=plt.cm.gray_r)
    nax3.get_xaxis().set_visible(False)
    nax3.get_yaxis().set_visible(False)

    color = plt.cm.rainbow(np.linspace(0, 1, n_colors))

    hrtree.t_n = 0
    def plotme(key):
        newSyn = np.zeros(shape)
        newSyn = core.synthesize_bubbles(newSyn, hrtree.clumps[key], hrtree.kernel, hrtree.noise, hrtree.delta)
        ct = nax.contour(newSyn.sum(axis=(0)), levels=[0.0], alpha=1.0, colors=[color[hrtree.t_n]])
        if print_numbers:
            plt.clabel(ct, fontsize=10, inline=1, fmt={0.0: str(key)})
        ct = nax2.contour(newSyn.sum(axis=(1)), levels=[0.0], alpha=1.0, colors=[color[hrtree.t_n]])
        if print_numbers:
            plt.clabel(ct, fontsize=10, inline=1, fmt={0.0: str(key)})
        ct = nax3.contour(newSyn.sum(axis=(2)).T, levels=[0.0], alpha=1.0, colors=[color[hrtree.t_n]])
        if print_numbers:
            plt.clabel(ct, fontsize=10, inline=1, fmt={0.0: str(key)})
        hrtree.t_n += 1

    def recu(tree):
        for (key, val) in tree.items():
            if hrtree.enabled[key]:
                if hrtree.display[key]:
                   plotme(key)
                if isinstance(val, dict):
                    recu(val)
    if hrtree.enabled[node]:
        if hrtree.display[node]:
            plotme(node)
    if isinstance(parent[node], dict):
        recu(parent[node])

    #i = 0
    #for b in bco:
    #    imask = (labels == b)
    #    npos = sol[imask]
    #    newSyn = np.zeros(shape)
    #    newSyn = acalib.core.synthesize_bubbles(newSyn, npos, mould, noise, delta)
    #    nax.contour(newSyn.sum(axis=(0)), levels=[0.0], alpha=1.0, colors=[color[i]])
    #    nax2.contour(newSyn.sum(axis=(1)), levels=[0.0], alpha=1.0, colors=[color[i]])
    #    nax3.contour(newSyn.sum(axis=(2)).T, levels=[0.0], alpha=1.0, colors=[color[i]])
    #    i += 1
    #plt.tight_layout()
    plt.show()
    #plt.figure()

    # nax4=fig.add_subplot(1,3,3)
    #nax4 = plt.subplot(111)
    #i = 0
    #for b in bco:
    #    imask = (labels == b)
    #    npos = sol[imask]
    #    newSyn = np.zeros(shape)
    #    newSyn = acalib.core.synthesize_bubbles(newSyn, npos, mould, noise, delta)
    #    nax.contour(newSyn.sum(axis=(0)), levels=[0.0], alpha=1.0, colors=[color[i]])
    #    nax2.contour(newSyn.sum(axis=(1)), levels=[0.0], alpha=1.0, colors=[color[i]])
    #    nax3.contour(newSyn.sum(axis=(2)).T, levels=[0.0], alpha=1.0, colors=[color[i]])
    #    nax4.plot(newSyn.sum(axis=(1, 2)), color=color[i])
    #    i += 1
    #plt.tight_layout()
    #plt.show()
    #return (bco)

@support_nddata
def visualize_image(data,wcs=None,unit=None,contour=False,cmap=None):
    """
    Plot 2D astronomical data.

    Parameters
    ------------
    data : numpy.ndarray or astropy.nddata.NDData or astropy.nddata.NDDataRef
        Astronomical image

    wcs : astropy.wcs.WCS
        World Coordinate System from the image (not needed if contained in NDData)

    unit : astropy.unit
        Image units (not needed if contained in NDData)

    contour : numpy.ndarray
        For plotting Contourns
    """
    if cmap is None:
        cmap=plt.cm.gray_r
    if wcs is None:
        if data.ndim != 2:
            log.info("Cannot visualize image data with no WCS and dimension != 2")
            return
        plt.imshow(data, origin='lower', cmap=cmap)
        cb = plt.colorbar()
        cb.ax.set_ylabel(unit)
    else:
        # Only work for spectral cubes... what else?
        if data.ndim == 3:
            try:
                freq_axis = axes._get_axis(wcs, "FREQ")
            except ValueError:
                log.warning("Data does not have a FREQ dimension to stack, aborting visualization")
                return
            wcs = wcs.dropaxis(wcs.naxis - freq_axis - 1)
            data = np.nansum(data, axis=(freq_axis))
        gax = plt.subplot(111, projection=wcs)
        _draw_image(data, gax, wcs, unit,cmap=cmap)
        g0 = gax.coords[0]
        g1 = gax.coords[1]
        g0.grid(color='yellow', alpha=0.5, linestyle='solid')
        g1.grid(color='yellow', alpha=0.5, linestyle='solid')
        if contour:
            arms=rms(data)
            dmax=data.max()
            crs=np.arange(1,dmax/arms)
            plt.contour(data,levels=arms*crs,alpha=0.5)
    plt.show()

import ipyvolume.pylab as ipvlab

@support_nddata
def visualize_volume(data,wcs=None,unit=None):
    if wcs is None:
        log.error("WCS is needed by this function")
        return

    if unit is None:
        log.error("Unit is needed by this function")
        return

    ipvlab.clear()

    labels = ["{} [{}]".format(axe, str(unit)) for axe, unit in
              zip(upi.axes_names(data, wcs), upi.axes_units(data, wcs))]

    ipvlab.xyzlabel(*labels[::-1])

    extent = upi.extent(data, wcs)
    minlim = extent[0]
    maxlim = extent[1]

    ipvlab.xlim(minlim[2].value, maxlim[2].value)
    ipvlab.ylim(minlim[1].value, maxlim[1].value)
    ipvlab.zlim(minlim[0].value, maxlim[0].value)

    ipvlab.style.use('dark')

    tf = ipvlab.transfer_function(level=[0.39, 0.54, 0.60], opacity=[0.2, 0.2, 0.2])
    vol = ipvlab.volshow(data, tf=tf, controls=False, level_width=0.1)

    ipvlab.gcf().width = 1024
    ipvlab.gcf().height = 456

    ipvlab.show()

    
@support_nddata
def visualize_contour3D(data,wcs=None,unit=None):
    pass
     # if wcs is None:
     #    log.error("WCS is needed by this function")
     # figure = mlab.figure('Contour Plot')
     # mesh=get_mesh(data)
     # xi,yi,zi=mesh
     # ranges=axes_ranges(data,wcs)
     # mmin = data.min()
     # mmax = data.max()
     # mlab.contour3d(xi,yi,zi,data,transparent=True,contours=10,opacity=0.5)
     # ax=mlab.axes(xlabel="VEL [km/s] ",ylabel="DEC [deg]",zlabel="RA [deg]",ranges=ranges,nb_labels=5)
     # ax.axes.label_format='%.3f'
     # mlab.colorbar(title='flux', orientation='vertical', nb_labels=5)
     # mlab.show()

# Specific plotting function
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

def visualize_rgb(rdata,gdata,bdata):
    wcs=rdata.wcs
    gax = plt.subplot(111, projection=wcs)
    sh = rdata.data.shape
    data = np.zeros((sh[0],sh[1],3))
    #plt.imshow(gdata.data)
    #plt.show()
    (rdata,_,_)=flux.standarize(rdata)
    (gdata,_,_)=flux.standarize(gdata)
    (bdata,_,_)=flux.standarize(bdata)
    data[:,:,0] = rdata.data/np.nanmax(rdata.data)
    data[:,:,1] = gdata.data/np.nanmax(gdata.data)
    data[:,:,2] = bdata.data/np.nanmax(bdata.data)
    plt.imshow(data, origin='lower')
    g0 = gax.coords[0]
    g1 = gax.coords[1]
    g0.set_axislabel(wcs.axis_type_names[0])
    g1.set_axislabel(wcs.axis_type_names[1])
    #cb = plt.colorbar()
    #cb.ax.set_ylabel(unit)
    g0 = gax.coords[0]
    g1 = gax.coords[1]
    g0.grid(color='yellow', alpha=0.5, linestyle='solid')
    g1.grid(color='yellow', alpha=0.5, linestyle='solid')
    plt.show()

from matplotlib import gridspec

# TODO: vmax can be computed by default
def show_image_grid(img_list,side,vmax,cmap=None):
    if cmap == None:
        cmap="gray_r"
    origin="lower"

    plt.figure(figsize=(10,10))
    gs = gridspec.GridSpec(side, side,wspace=0.0, hspace=0.0)
    for i in range(side*side):
        ax  = plt.subplot(gs[i])
        ax.imshow(img_list[i],origin=origin,vmax=vmax,vmin=0,cmap=cmap)
        ax.get_xaxis().set_visible(False)
        ax.get_yaxis().set_visible(False)
    plt.show()

