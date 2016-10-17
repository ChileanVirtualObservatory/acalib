import numpy as np

def get_mesh(data,lower=None,upper=None):
    """ Create an index mesh """
    sl=slab(data,lower,upper)
    dim=data.ndim
    slices=[]
    for i in range(dim):
       slices.append(slice(sl[i].start,sl[i].stop))
    retval=np.mgrid[slices]
    return retval

def to_features(data,lower=None,upper=None):
    msh=get_mesh(data,lower,upper)
    dim=data.ndim
    ii=np.empty((dim,int(msh.size/dim)))
    for i in range(dim):
       ii[dim-i-1]=msh[i].ravel()
    return ii

support_nddata
def get_velocities(data,wcs=None,fqi=None,restfrq=None):
    if wcs is None:
        log.error("A world coordinate system (WCS) is needed")
        return None
    if fqi is None:
        return None
    if restfrq is None:
        restfrq=wcs.wcs.restfrq*u.Hz
    dim=wcs.wcs.spec
    idx=np.zeros((fqi.size,data.ndim))
    idx[:,dim]=fqi
    vals=wcs.all_pix2world(idx,0)
    eq=u.doppler_radio(restfrq)
    vec=vals[:,dim]*u.Hz
    return vec.to(u.km/u.s, equivalencies=eq)

# TODO: extend to n-dimensions (only works for 3)
@support_nddata
def axes_ranges(data,wcs=None,lower=None,upper=None):
    """ Get axes extent (transforms freq to velocity!) """
    if lower==None:
        lower=[0,0,0]
    if upper==None:
        upper=data.shape
    lower=lower[::-1]
    lwcs=wcs.wcs_pix2world([lower], 0)
    lwcs=lwcs[0]
    upper=upper[::-1]
    uwcs=wcs.wcs_pix2world([upper], 0)
    uwcs=uwcs[0]
    lfreq=lwcs[2]*u.Hz
    ufreq=uwcs[2]*u.Hz
    rfreq=wcs.wcs.restfrq*u.Hz
    eq= u.doppler_radio(rfreq)
    lvel=lfreq.to(u.km/u.s, equivalencies=eq)
    uvel=ufreq.to(u.km/u.s, equivalencies=eq)
    ranges=[lvel.value,uvel.value,lwcs[1],uwcs[1],lwcs[0],uwcs[0]]
    return ranges

#TODO: try to merge with axes_ranges and get_velocities!
@support_nddata
def axis_range(data,wcs,axis):
    lower=wcs.wcs_pix2world([[0,0,0]], 0) - wcs.wcs.cdelt/2.0
    shape=data.shape
    shape=[shape[::-1]]
    upper=wcs.wcs_pix2world(shape, 1) + wcs.wcs.cdelt/2.0
    return (lower[0][axis],upper[0][axis])



### DEPRECATED ####
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

