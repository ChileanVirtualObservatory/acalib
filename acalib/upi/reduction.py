from astropy.nddata import support_nddata, NDDataRef
from astropy.table import Table
from astropy import log
import numpy as np
import astropy.units as u
from acalib import core
from acalib.upi.axes import spectral_velocities



@support_nddata
def _moment(data,order,wcs=None,mask=None,unit=None,restfrq=None):
    if wcs is None:
        log.error("A world coordinate system (WCS) is needed")
        return None
    data=core.fix_mask(data,mask)
    dim=wcs.wcs.spec
    rdim=data.ndim - 1 - dim
    v=spectral_velocities(data,wcs,fqis=np.arange(data.shape[rdim]),restfrq=restfrq)
    v=v.value
    m0=data.sum(axis=rdim)
    if order==0:
        mywcs=wcs.dropaxis(dim)
        return NDDataRef(m0.data, uncertainty=None, mask=m0.mask,wcs=mywcs, meta=None, unit=unit)
    #mu,alpha=np.average(data,axis=rdim,weights=v,returned=True)
    mu,alpha=np.ma.average(data,axis=rdim,weights=v,returned=True)
    m1=alpha*mu/m0
    if order==1:
        mywcs=wcs.dropaxis(dim)
        return NDDataRef(m1.data, uncertainty=None, mask=m1.mask,wcs=mywcs, meta=None, unit=u.km/u.s)
    v2=v*v
    var,beta=np.ma.average(data,axis=rdim,weights=v2,returned=True)
    #var,beta=data.average(axis=rdim,weights=v2,returned=True)
    m2=np.sqrt(beta*var/m0 - m1*m1)
    if order==2:
        mywcs=wcs.dropaxis(dim)
        return NDDataRef(m2.data, uncertainty=None, mask=m2.mask,wcs=mywcs, meta=None, unit=u.km*u.km/u.s/u.s)
    log.error("Order not supported")
    return None

# Should return a NDData
@support_nddata
def moment0(data,wcs=None,mask=None,unit=None,restfrq=None):
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
    return _moment(data,0,wcs,mask,unit,restfrq)

@support_nddata
def moment1(data,wcs=None,mask=None,unit=None,restfrq=None):
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
    return _moment(data,1,wcs,mask,unit,restfrq)

@support_nddata
def moment2(data,wcs=None,mask=None,unit=None,restfrq=None):
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
    return _moment(data,2,wcs,mask,unit,restfrq)


# TODO: Fix this function, is not working correctly
@support_nddata
def spectra(data,wcs=None,mask=None,unit=None,restrict=None):
    """


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
        restrict : boolean


        Returns
        -------
        result: astropy.nddata.NDData
            Moment 2 of the data cube

    """
    if restrict is None:
        #Create NDD and WCS change...
        return core.integrate(data,axis=(1,2))
    else:
        log.error("Not Implemented Yet!")

        # Get 1 pixel aperture
        aperture=np.abs(wcs.celestial.wcs.cdelt[0])*u.deg
    #if position.unit == u.pix and aperture.unit == u.pix:
    #    # TODO:  Here is the nasty part
    #    lb=np.array([0,            position[1].value - aperture.value, position[0].value - aperture.value])
    #    ub=np.array([data.shape[2],position[1].value + aperture.value, position[0].value + aperture.value])
    #else:
    #    log.error("Not Implemented Yet!")
    #specview=data[slab(data,lb,ub)]
    #return specview.sum(axis=(1,2))
