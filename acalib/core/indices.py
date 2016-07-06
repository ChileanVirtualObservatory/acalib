import numpy as np


def fix_limits(data,vect):
    """ Fix vect index to be inside data """ 
    if isinstance(vect,tuple):
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

def world_window_to_index(data,wcs,center,window):
    ld=np.rint(wcs.wcs_world2pix([center-window],0))
    lu=np.rint(wcs.wcs_world2pix([center+window],0))
    lower=np.array([ld,lu]).min(axis=0)
    upper=np.array([ld,lu]).max(axis=0)
    lower=fix_limits(data,lower[0][::-1])
    upper=fix_limits(data,upper[0][::-1])
    return (lower,upper)


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
       ii[dim-i-1]=msh[i].ravel()
       ii[dim-i-1]=msh[i].ravel()
    return ii


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

