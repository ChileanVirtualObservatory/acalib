import numpy as np
import astropy.units as u
from astropy.nddata import support_nddata, NDData
from astropy.table import Table, Column
from astropy import log

import core
from upi.formatting import _unitize, _world_table_creator

@support_nddata
def cut(data, wcs=None, mask=None, unit=None, lower=None, upper=None):
    # Check for NDDataSlicing... maybe this is already done by astropy.nddata package.
    mslab = core.slab(data, lower, upper)
    scube = data[mslab]
    newwcs = wcs.slice(mslab, numpy_order=True)
    return NDData(scube, wcs=newwcs, unit=unit)

@support_nddata
def axes_names(data,wcs=None):
    if wcs is None:
        log.error("A world coordinate system (WCS) is needed")
        return None
    return np.array(wcs.axis_type_names)[::-1]


@support_nddata
def extent(data,wcs=None,lower=None,upper=None):
    """ Get axes extent  """
    #TODO: These can be a decorator
    if wcs is None:
        log.error("A world coordinate system (WCS) is needed")
        return None
    if lower==None:
        lower=np.zeros(data.ndim)
    if upper==None:
        upper=data.shape
    idx=[lower,upper]
    idx_f  = np.fliplr(idx)
    values = wcs.wcs_pix2world(idx_f, 0)
    values = np.fliplr(values)
    return (_unitize(values[0],wcs),_unitize(values[1],wcs))

@support_nddata
def center(data,wcs=None):
    """ Get center of the data"""
    #TODO: These can be a decorator
    if wcs is None:
        log.error("A world coordinate system (WCS) is needed")
        return None
    val=wcs.wcs.crval[::-1]
    return _unitize(val,wcs)

@support_nddata
def axes_units(data,wcs=None):
    """ Get units of the axes"""
    #TODO: These can be a decorator (
    if wcs is None:
        log.error("A world coordinate system (WCS) is needed")
        return None
    uvec=np.array(wcs.wcs.cunit)[::-1]
    return uvec

@support_nddata
def resolution(data,wcs=None):
    """ Get the resolution of data"""
    if wcs is None:
        log.error("A world coordinate system (WCS) is needed")
        return None
    val=wcs.wcs.cdelt[::-1]
    return _unitize(val,wcs)

@support_nddata
def spectral_velocities(data,wcs=None,fqs=None,fqis=None,restfrq=None):
    """ Get the spectral velocities from frequencies fqs given a rest frequency (by default search for it in the WCS). If fqs is None, then frequencies indices (fqis) need to be given.
    """
    if wcs is None:
        log.error("A world coordinate system (WCS) is needed")
        return None
    if restfrq is None:
        restfrq=wcs.wcs.restfrq*u.Hz
    if fqs is None:
        if fqis is None:
           return None
        dim=wcs.wcs.spec
        idx=np.zeros((fqis.size,data.ndim))
        idx[:,dim]=fqis
        vals=wcs.all_pix2world(idx,0)
        fqs=vals[:,dim]*u.Hz
    eq=u.doppler_radio(restfrq)
    return fqs.to(u.km/u.s, equivalencies=eq)


@support_nddata
def features(data,wcs=None,lower=None,upper=None):
    """ Creates an array with WCS axea in features format """
    if wcs is None:
        log.error("A world coordinate system (WCS) is needed")
        return None
    ii=core.index_features(data,lower,upper)
    f=wcs.wcs_pix2world(ii.T,0)
    return _world_table_creator(f,wcs)


@support_nddata
# TODO: Consider using "box" structure rather than up and low
def opening(data,center,window,wcs=None):
    """ Field of view (center +- window) converted to indices """
    if wcs is None:
        log.error("A world coordinate system (WCS) is needed")
        return None
    off_low=center-window
    off_low = np.array([x.value for x in off_low])
    off_up = center+window
    off_up = np.array([x.value for x in off_up])
    #dim = len(center.colnames)
    #off_low = np.array([center[0][i] - window[0][i] for i in range(dim)])
    #off_up  = np.array([center[0][i] + window[0][i] for i in range(dim)])
    ld=np.rint(wcs.wcs_world2pix([off_low[::-1]],0))
    lu=np.rint(wcs.wcs_world2pix([off_up[::-1]],0))
    lower=np.array([ld,lu]).min(axis=0)
    upper=np.array([ld,lu]).max(axis=0)
    lower=core.fix_limits(data,lower[0][::-1])
    upper=core.fix_limits(data,upper[0][::-1])
    return (lower,upper)
    #values=np.vstack((lower,upper))
    #return _pix_table_creator(values,wcs)
    


### DEPRECATED ####
#TODO: try to merge with axes_ranges and get_velocities!
#@support_nddata
#def axis_range(data,wcs,axis):
#    lower=wcs.wcs_pix2world([[0,0,0]], 0) - wcs.wcs.cdelt/2.0
#    shape=data.shape
#    shape=[shape[::-1]]
#    upper=wcs.wcs_pix2world(shape, 1) + wcs.wcs.cdelt/2.0
#    return (lower[0][axis],upper[0][axis])
#def ndslice(ndd, lower, upper):
#    """ 
#    N-Dimensional slicing.
#    
#    Arguments:
#        ndd   -- an astropy.nddata.NDDataArray object.
#        lower -- n-dimensional point as an n-tuple.
#        upper -- n-dimensional point as an n-tuple.
#    
#    Returns:
#        A sliced astropy.nddata.NDDataArray object.
#        
#    """
#    lower = lower if lower is not None else np.zeros(ndd.ndim)
#    upper = upper if upper is not None else ndd.shape
#    return ndd[[slice(min(a,b), max(a,b)+1) for a,b in zip(lower, upper)]]
#
#def adjust_index(relative, origin):
#    """
#    Adjusts an index relative to a subarray to an absolute
#    index in the superarray.
#    
#    Arguments:
#        origin   -- an n-dimensional index of a point as an n-tuple.
#                    It should be the origin from which the relative
#                    index was computed.
#        relative -- an n-dimensional index of a point as an n-tuple.
#                    The index to be adjusted.
#    
#    Returns:
#        The relative index adjusted to the superarray as an n-tuple.
#    """
#    return tuple(np.array(origin) + np.array(relative))

#def index_of_max(ndd, lower=None, upper=None):
#    """ 
#    Index of maximum value in an m-dimensional subarray from 
#    an n-dimensional array, specified by lower and upper.
#    
#    Arguments:
#        ndd   -- an astropy.nddata.NDDataArray object.
#        lower -- n-dimensional point as an n-tuple.
#        upper -- n-dimensional point as an n-tuple.
#    
#    Returns:
#        A tuple with the maximum value found in the m-dimensional
#        subarray and its index in the n-dimensional superarray.
#        
#    """
#    ndd = ndslice(ndd, lower, upper)
#    index = np.unravel_index(ndd.data.argmax(), ndd.data.shape)
#    value = ndd.data[index]
#    return (value, adjust_index(index, lower))

#def index_of_min(ndd, lower=None, upper=None):
#    """ 
#    Index of minimum value in an m-dimensional subarray from 
#    an n-dimensional array, specified by lower and upper.
#    
#    Arguments:
#        ndd   -- an astropy.nddata.NDDataArray object.
#        lower -- n-dimensional point as an n-tuple.
#        upper -- n-dimensional point as an n-tuple.
#    
#    Returns:
#        A tuple with the minimum value found in the m-dimensional
#        subarray and its index in the n-dimensional superarray.
#        
#    """
#    ndd = ndslice(ndd, lower, upper)
#    index = np.unravel_index(ndd.data.argmin(), ndd.data.shape)
#    value = ndd.data[index]
#    return (value, adjust_index(index, lower))

