import numpy as np
import astropy.units as u
from astropy.nddata import support_nddata, NDData

def fix_limits(data,vect):
    """ Fix vect index to be inside data """
    if isinstance(vect,(tuple,list)):
        vect=np.array(vect)
    vect=vect.astype(int)
    low=vect < 0
    up=vect > data.shape
    if vect.any():
        vect[low]=0
    if vect.any():
        vect[up]=np.array(data.shape)[up]
    return vect


def slab(data,lower=None,upper=None):
    """ Obtain the n-dimensional slab from lower to upper (i.e. slab is a vector of slices)"""
    if lower is None:
        lower=np.zeros(data.ndim)
    if upper is None:
        upper=data.shape
    lower=fix_limits(data,lower)
    upper=fix_limits(data,upper)
    m_slab=[]
    for i in range(data.ndim):
       m_slab.append(slice(lower[i],upper[i]))
    return m_slab


def matching_slabs(data,flux,lower,upper):
    """ Obtain the matching data and flux slabs from lower to upper while fixing the limits"""
    data_slab=slab(data,lower,upper)
    flow=np.zeros(flux.ndim)
    fup=np.array(flux.shape)
    for i in range(data.ndim):
       if data_slab[i].start == 0:
          flow[i] = flux.shape[i] - data_slab[i].stop
       if data_slab[i].stop == data.shape[i]:
          fup[i] = data_slab[i].stop - data_slab[i].start
    flux_slab=slab(flux,flow,fup)
    return data_slab,flux_slab


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
def extent(data,wcs=None,lower=None,upper=None):
    """ Get axes extent  """
    if wcs is None:
        log.error("A world coordinate system (WCS) is needed")
        return None
    if lower==None:
        lower=np.zeros(data.ndim)
    if upper==None:
        upper=data.shape
    uvec=np.array(wcs.wcs.cunit)
    lower=lower[::-1]
    upper=upper[::-1]
    lwcs=wcs.wcs_pix2world([lower], 0)
    lwcs=lwcs[0]*uvec
    uwcs=wcs.wcs_pix2world([upper], 0)
    uwcs=uwcs[0]*uvec
    lwcs=lwcs[::-1]
    uwcs=uwcs[::-1]
    ranges=np.array([lwcs,uwcs]).T.ravel()
    return ranges

@support_nddata
def center(data,wcs=None):
    """ Get center of the data"""
    if wcs is None:
        log.error("A world coordinate system (WCS) is needed")
        return None
    uvec=np.array(wcs.wcs.cunit)
    val=np.array(wcs.wcs.crval)*uvec
    return val[::-1]

@support_nddata
def resolution(data,wcs=None):
    """ Get the resolution of data"""
    if wcs is None:
        log.error("A world coordinate system (WCS) is needed")
        return None
    uvec=np.array(wcs.wcs.cunit)
    val=np.array(wcs.wcs.cdelt)*uvec
    return val[::-1]

@support_nddata
def axes_names(data,wcs=None):
    if wcs is None:
        log.error("A world coordinate system (WCS) is needed")
        return None
    return np.array(wcs.axis_type_names)[::-1]

@support_nddata
def index_mesh(data,lower=None,upper=None):
    """ Create an meshgrid from indices """
    sl=slab(data,lower,upper)
    dim=data.ndim
    slices=[]
    for i in range(dim):
       slices.append(slice(sl[i].start,sl[i].stop))
    retval=np.mgrid[slices]
    return retval

@support_nddata
def index_features(data,lower=None,upper=None):
    """ Creates an array with indices in features format """
    msh=index_mesh(data,lower,upper)
    dim=data.ndim
    ii=np.empty((dim,int(msh.size/dim)))
    for i in range(dim):
       ii[dim-i-1]=msh[i].ravel()
    return ii


@support_nddata
def features(data,wcs=None,lower=None,upper=None):
    """ Creates an array with WCS axea in features format """
    if wcs is None:
        log.error("A world coordinate system (WCS) is needed")
        return None
    ii=index_features(data,lower,upper)
    f=wcs.wcs_pix2world(ii.T,0)
    ff=np.multiply(f,np.array(wcs.wcs.cunit))
    return np.fliplr(ff)


@support_nddata
# TODO: Consider using "box" structure rather than up and low
def opening(data,center,window,wcs=None):
    """ Field of view (center +- window) converted to indices """
    if wcs is None:
        log.error("A world coordinate system (WCS) is needed")
        return None
    off_low=center-window
    off_up = center+window
    off_low = np.array([x.value for x in off_low])
    off_up = np.array([x.value for x in off_up])
    ld=np.rint(wcs.wcs_world2pix([off_low[::-1]],0))
    lu=np.rint(wcs.wcs_world2pix([off_up[::-1]],0))
    lower=np.array([ld,lu]).min(axis=0)
    upper=np.array([ld,lu]).max(axis=0)
    lower=fix_limits(data,lower[0][::-1])
    upper=fix_limits(data,upper[0][::-1])
    return (lower,upper)


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

