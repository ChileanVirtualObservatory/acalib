import numpy as np
from utils import ndslice, adjust_index

def index_of_max(ndd, lower=None, upper=None):
    """ 
    Index of maximum value in an m-dimensional subarray from 
    an n-dimensional array, specified by lower and upper.
    
    Arguments:
        ndd   -- an astropy.nddata.NDDataArray object.
        lower -- n-dimensional point as an n-tuple.
        upper -- n-dimensional point as an n-tuple.
    
    Returns:
        A tuple with the maximum value found in the m-dimensional
        subarray and its index in the n-dimensional superarray.
        
    """
    ndd = ndslice(ndd, lower, upper)
    index = np.unravel_index(ndd.data.argmax(), ndd.data.shape)
    value = ndd.data[index]
    return (value, adjust_index(index, lower))

def index_of_min(ndd, lower=None, upper=None):
    """ 
    Index of minimum value in an m-dimensional subarray from 
    an n-dimensional array, specified by lower and upper.
    
    Arguments:
        ndd   -- an astropy.nddata.NDDataArray object.
        lower -- n-dimensional point as an n-tuple.
        upper -- n-dimensional point as an n-tuple.
    
    Returns:
        A tuple with the minimum value found in the m-dimensional
        subarray and its index in the n-dimensional superarray.
        
    """
    ndd = ndslice(ndd, lower, upper)
    index = np.unravel_index(ndd.data.argmin(), ndd.data.shape)
    value = ndd.data[index]
    return (value, adjust_index(index, lower))

def index_to_world(ndd, index):
    pass

def index_as_features(ndd, lower, upper):
    pass

def index_from_window(ndd, wcs_center, wcs_window):
    pass

def index_fix(ndd, index):
    pass

def index_mesh(ndd):
    pass
