from astropy.nddata import support_nddata, NDData
from astropy.table import Table
from astropy import log
import numpy as np
import astropy.units as u
from acalib import core
from acalib.upi.axes import spectral_velocities



# TODO: pack/unpack with the NDData decorator?
def _moment(data,order,wcs=None,mask=None,unit=None,restfrq=None): 
    # TODO: Decorator?
    if wcs is None:
        log.error("A wo rld coordinate system (WCS) is needed")
        return None 
    data=core.fix_mask(data,mask)
    dim=wcs.wcs.spec 
    rdim=data.ndim - 1 - dim 
    v=spectral_velocities(data,wcs,fqis=np.arange(data.shape[rdim]),restfrq=restfrq) 
    v=v.value
    m0=data.sum(axis=rdim) 
    if order==0: 
        mywcs=wcs.dropaxis(dim) 
        return NDData(m0.data, uncertainty=None, mask=m0.mask,wcs=mywcs, meta=None, unit=unit) 
    #mu,alpha=np.average(data,axis=rdim,weights=v,returned=True) 
    mu,alpha=np.ma.average(data,axis=rdim,weights=v,returned=True) 
    m1=alpha*mu/m0 
    if order==1: 
        mywcs=wcs.dropaxis(dim) 
        return NDData(m1.data, uncertainty=None, mask=m1.mask,wcs=mywcs, meta=None, unit=u.km/u.s) 
    v2=v*v 
    var,beta=np.ma.average(data,axis=rdim,weights=v2,returned=True) 
    #var,beta=data.average(axis=rdim,weights=v2,returned=True) 
    m2=np.sqrt(beta*var/m0 - m1*m1) 
    if order==2: 
        mywcs=wcs.dropaxis(dim) 
        return NDData(m2.data, uncertainty=None, mask=m2.mask,wcs=mywcs, meta=None, unit=u.km*u.km/u.s/u.s) 
    log.error("Order not supported") 
    return None
        
# Should return a NDData 
@support_nddata 
def moment0(data,wcs=None,mask=None,unit=None,restfrq=None): 
    return _moment(data,0,wcs,mask,unit,restfrq) 
 
@support_nddata 
def moment1(data,wcs=None,mask=None,unit=None,restfrq=None): 
    return _moment(data,1,wcs,mask,unit,restfrq) 
 
@support_nddata 
def moment2(data,wcs=None,mask=None,unit=None,restfrq=None): 
    return _moment(data,2,wcs,mask,unit,restfrq) 


# TODO: Fix this function, is not working...
@support_nddata
def spectra(data,wcs=None,mask=None,unit=None,position=None,aperture=None):
    if position is None:
        # Get celestial center
        position=wcs.celestial.wcs.crval*u.deg
    if aperture is None:
        # Get 1 pixel aperture
        aperture=np.abs(wcs.celestial.wcs.cdelt[0])*u.deg
    if position.unit == u.pix and aperture.unit == u.pix:
        # TODO:  Here is the nasty part
        lb=np.array([0,            position[1].value - aperture.value, position[0].value - aperture.value])
        ub=np.array([data.shape[2],position[1].value + aperture.value, position[0].value + aperture.value])
    else:
        log.error("Not Implemented Yet!")
    specview=data[slab(data,lb,ub)]
    return specview.sum(axis=(1,2))
