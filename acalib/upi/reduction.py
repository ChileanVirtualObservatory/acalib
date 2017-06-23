from astropy.nddata import support_nddata
from astropy.table import Table
from astropy import log
import numpy as np
import astropy.units as u

from ipywidgets import interact

from acalib import core
from acalib.io import graph
from acalib.upi import axes
from acalib.upi.data import Data
from acalib.upi.axes import spectral_velocities


@support_nddata
def _moment(data, order, wcs=None, mask=None, unit=None, meta=None, restfrq=None):
    if wcs is None:
        log.error("A world coordinate system (WCS) is needed")
        return None
    data = core.fix_mask(data, mask)
    dim = wcs.wcs.spec
    rdim = data.ndim - 1 - dim
    m0 = data.sum(axis=rdim)

    v = spectral_velocities(data, wcs, fqis=np.arange(data.shape[rdim]), restfrq=restfrq)
    v = v.value
    # mu,alpha=np.average(data,axis=rdim,weights=v,returned=True)
    if order == 0:
        mywcs = wcs.dropaxis(dim)
        return Data(m0.data, uncertainty=None, mask=m0.mask, wcs=mywcs, meta=meta, unit=unit)
    mu, alpha = np.ma.average(data, axis=rdim, weights=v, returned=True)
    m1 = alpha * mu / m0
    if order == 1:
        mywcs = wcs.dropaxis(dim)
        return Data(m1.data, uncertainty=None, mask=m1.mask, wcs=mywcs, meta=meta, unit=u.km / u.s)
    v2 = v * v
    var, beta = np.ma.average(data, axis=rdim, weights=v2, returned=True)
    # var,beta=data.average(axis=rdim,weights=v2,returned=True)
    m2 = np.sqrt(beta * var / m0 - m1 * m1)
    if order == 2:
        mywcs = wcs.dropaxis(dim)
        return Data(m2.data, uncertainty=None, mask=m2.mask, wcs=mywcs, meta=meta, unit=u.km * u.km / u.s / u.s)
    log.error("Order not supported")
    return None

@support_nddata
def moment0(data, wcs=None, mask=None, unit=None, meta=None, restfrq=None):
    """
        Calculate moment 0 from a data cube.

        Parameters
        ----------
        data : (M,N,Z) numpy.ndarray or astropy.nddata.NDData or astropy.nddata.NDDataRef
            Astronomical data cube.
        wcs : astropy.wcs.wcs.WCS
            World Coordinate System to use.
        mask : numpy.ndarray
            Mask for data.
        unit : astropy.units.Unit
            Astropy unit (http://docs.astropy.org/en/stable/units/).
        restfrq : astropy.units.quantity.Quantity
            Rest frequency

        Returns
        -------
        result: astropy.nddata.NDDataRef
            Moment 0 of the data cube

    """
    return _moment(data, 0, wcs, mask, unit, meta, restfrq)


@support_nddata
def moment1(data, wcs=None, mask=None, unit=None, meta=None, restfrq=None):
    """
        Calculate moment 1 from a data cube.

        Parameters
        ----------
        data : (M,N,Z) numpy.ndarray or astropy.nddata.NDData
            Astronomical data cube.
        wcs : astropy.wcs.wcs.WCS
            World Coordinate System to use
        mask : numpy.ndarray
            Mask for data.
        unit : astropy.units.Unit
            Astropy unit (http://docs.astropy.org/en/stable/units/)
        restfrq : astropy.units.quantity.Quantity
            Rest frequency

        Returns
        -------
        result: astropy.nddata.NDData
            Moment 1 of the data cube

    """
    return _moment(data, 1, wcs, mask, unit, meta, restfrq)


@support_nddata
def moment2(data, wcs=None, mask=None, unit=None, meta=None, restfrq=None):
    """
        Calculate moment 2 from a data cube.

        Parameters
        ----------
        data : (M,N,Z) numpy.ndarray or astropy.nddata.NDData or astropy.nddata.NDDataRef
            Astronomical data cube.
        wcs : astropy.wcs.wcs.WCS
            World Coordinate System to use
        mask : numpy.ndarray
            Mask for data.
        unit : astropy.units.Unit
            Astropy unit (http://docs.astropy.org/en/stable/units/)
        restfrq : astropy.units.quantity.Quantity
            Rest frequency

        Returns
        -------
        result: astropy.nddata.NDDataRef
            Moment 2 of the data cube

    """
    return _moment(data, 2, wcs, mask, unit, meta, restfrq)


@support_nddata
def select_region(data, wcs=None, mask=None, unit=None, meta=None,ra_1=None,dec_1=None,ra_2=None,dec_2=None,interactive=False):
    ra_idx = axes._get_axis(wcs,'RA')
    dec_idx = axes._get_axis(wcs,'DEC')
    ext1,ext2 = axes.extent(data, wcs=wcs)
    if ra_1 is None:
        ra_1 = ext1[ra_idx]
        dec_1 = ext1[dec_idx]
    if ra_2 is None:
        ra_2 = ext2[ra_idx]
    if dec_2 is None:
        dec_2 = ext2[dec_idx]
    # All in arcsecs

    ra_1=  ra_1.to("arcsec").value
    ra_2 = ra_2.to("arcsec").value
    dec_1 = dec_1.to("arcsec").value
    dec_2 = dec_2.to("arcsec").value
    min_ra = min(ra_1, ra_2)
    max_ra = max(ra_1, ra_2)
    min_dec = min(dec_1, dec_2)
    max_dec = max(dec_1, dec_2)

    def display_region(RA1=ra_1, RA2 = ra_2, DEC1=dec_1, DEC2=dec_2):
        if data.ndim == 3:
            lower=[0,0,0]
            upper=[0,0,0]
            freq_idx = axes._get_axis(wcs,'FREQ')
            lower[freq_idx] = ext1[freq_idx].value
            upper[freq_idx] = ext2[freq_idx].value
        else:
            lower=[0,0]
            upper=[0,0]
        lower[ra_idx] = (RA1*u.arcsec).to("deg").value
        upper[ra_idx] = (RA2*u.arcsec).to("deg").value
        lower[dec_idx] = (DEC1*u.arcsec).to("deg").value
        upper[dec_idx] = (DEC2*u.arcsec).to("deg").value
        upper = upper[::-1]
        lower = lower[::-1]
        ii=np.array([lower,upper])
        f = wcs.wcs_world2pix(ii, 0)
        values = np.fliplr(f)
        [lower,upper]=values.astype("int")
        if interactive == True:
            graph.show_subcube(data,wcs,meta,mask,unit,lower,upper)
        return np.array([lower,upper])

    if interactive == True:
        res=interact(display_region, RA1=(min_ra,max_ra), RA2=(min_ra,max_ra),
              DEC1=(min_dec,max_dec),DEC2=(min_dec,max_dec))
        result=res.widget.result
    else:
        result=display_region()
    return result

@support_nddata
def select_band(data, wcs=None, mask=None, unit=None, meta=None,freq_1=None,freq_2=None,interactive=False):
    if data.ndim != 3:
        log.warning("Bandwidth selection available only in 3D data ")
        return None
    ra_idx = axes._get_axis(wcs,'RA')
    dec_idx = axes._get_axis(wcs,'DEC')
    freq_idx = axes._get_axis(wcs,'FREQ')
    ext1,ext2 = axes.extent(data, wcs=wcs)
    if freq_1 is None:
        freq_1=ext1[freq_idx]
    if freq_2 is None:
        freq_2=ext2[freq_idx]
    freq_1=freq_1.value
    freq_2 = freq_2.value
    #ra_1 = ext1[ra_idx]
    #dec_1 = ext1[dec_idx]
    #ra_2 = ext2[ra_idx]
    #dec_2 = ext2[dec_idx]
    # All in arcsecs
    #ra_1=ra_1.value
    #ra_2 = ra_2.value
    #dec_1 = dec_1.value
    #dec_2 = dec_2.to("arcsec").value
    #min_ra = min(ra_1, ra_2)
    #max_ra = max(ra_1, ra_2)
    #min_dec = min(dec_1, dec_2)
    #max_dec = max(dec_1, dec_2)

    def display_region(FREQ1=freq_1,FREQ2=freq_2):
        lower=[0,0,0]
        lower[ra_idx]=ext1[ra_idx].value
        lower[dec_idx] = ext1[dec_idx].value
        lower[freq_idx] = FREQ1
        upper=[0,0,0]
        upper[ra_idx] = ext2[ra_idx].value
        upper[dec_idx] = ext2[dec_idx].value
        upper[freq_idx] = FREQ2
        upper = upper[::-1]
        lower = lower[::-1]
        ii=np.array([lower,upper])
        f = wcs.wcs_world2pix(ii, 0)
        values = np.fliplr(f)
        [lower,upper]=values.astype("int")
        if interactive == True:
            graph.show_subcube(data,wcs,meta,mask,unit,lower,upper)
        return np.array([lower,upper])

    if interactive == True:
        res=interact(display_region, FREQ1=(freq_1,freq_2),FREQ2=(freq_1,freq_2))
        result=res.widget.result
    else:
        result=display_region();
    return result


@support_nddata
def cut(data, wcs=None, mask=None, unit=None, meta=None, region=None):
    """
        Get a cut of the cube.

        Parameters
        ----------
        data : (M,N) or (M,N,Z) numpy.ndarray or astropy.nddata.NDData or astropy.nddata.NDDataRef
            Astronomical data cube.
        wcs : astropy.wcs.wcs.WCS
            World Coordinate System to use.
        mask : numpy.ndarray
            Mask for data.
        unit : astropy.units.Unit
            Astropy unit (http://docs.astropy.org/en/stable/units/).
        lower : tuple
            Start index from where to cut.
        upper : tuple
            Index to end cut.

        Returns
        -------
         result: acalib.upi.AData.
            Data cut from lower to upper.
    """
    # Check for NDDataSlicing... maybe this is already done by astropy.nddata package.
    if region is None:
        lower = None
        upper = None
    else:
        (lower, upper) = region
    mslab = core.slab(data, lower, upper)
    scube = data[mslab]
    if mask is None:
        smask=None
    else:
        smask = mask[mslab]
    newwcs = wcs.slice(mslab, numpy_order=True)
    return Data(scube, wcs=newwcs, unit=unit, mask=smask,meta=meta)

# TODO: Fix this function, is not working correctly (rezise the wcs... uff)
# TODO: Deprecated!
# @support_nddata
# def spectra(data, wcs=None, mask=None, unit=None, restrict=None):
#     """
#
#         Parameters
#         ----------
#         data : (M,N,Z) numpy.ndarray or astropy.nddata.NDData or astropy.nddata.NDDataRef
#             Astronomical data cube.
#         wcs : astropy.wcs.wcs.WCS
#             World Coordinate System to use
#         mask : numpy.ndarray
#             Mask for data.
#         unit : astropy.units.Unit
#             Astropy unit (http://docs.astropy.org/en/stable/units/)
#         restrict : boolean
#     """
#     if restrict is None:
#         # Create NDD and WCS change...
#         return core.integrate(data, axis=(1, 2))
#     else:
#         log.error("Not Implemented Yet!")
#
#         # Get 1 pixel aperture
#         aperture = np.abs(wcs.celestial.wcs.cdelt[0]) * u.deg
#         # if position.unit == u.pix and aperture.unit == u.pix:
#         #    # TODO:  Here is the nasty part
#         #    lb=np.array([0,            position[1].value - aperture.value, position[0].value - aperture.value])
#         #    ub=np.array([data.shape[2],position[1].value + aperture.value, position[0].value + aperture.value])
#         # else:
#         #    log.error("Not Implemented Yet!")
#         # specview=data[slab(data,lb,ub)]
#         # return specview.sum(axis=(1,2))
#
