import numpy as np

def fix_mask(data, mask):
    """

    Parameters
    ----------
    data : numpy.ndarray or numpy.ma.MaskedArray
        Astronomical data cube.
    mask : numpy.ndarray
        Boolean that will be applied.
    
    Returns
    -------
    result : numpy.ma.MaskedArray
        Masked astronomical data cube.
    """

    ismasked = isinstance(data, np.ma.MaskedArray)
    if ismasked and mask is None:
        return data
    else:
        return np.ma.MaskedArray(data, mask)

def fix_limits(data, vect):
    """ 
    Fix vect index to be inside data

    Parameters
    ----------
    data : numpy.ndarray or numpy.ma.MaskedArray
        Astronomical data cube.
    vect : tuple, list or numpy.ndarray
        Array with the indexes to be fixed.

    Returns
    -------
    result : numpy.ndarray
        Fixed array of indexes.
    """

    if isinstance(vect, (tuple, list)):
        vect = np.array(vect)
    vect = vect.astype(int)
    low = vect < 0
    up = vect > data.shape
    if vect.any():
        vect[low] = 0
    if vect.any():
        vect[up] = np.array(data.shape)[up]
    return vect


def slab(data, lower=None, upper=None):
    """
    Obtain the n-dimensional slab from lower to upper (i.e. slab is a vector of slices)
    
    Parameters
    ----------
    data : numpy.ndarray
        Atronomical data cube.
    lower : 3-tuple (default=None)
        Lower coordinates for the subcube.

    upper : 3-tuple (default=None)
        Upper coordinates for the subcube.


    Returns
    -------
    result : list
       list of slices using lower and upper coordinates to create a subcube.
    """

    if lower is None:
        lower = np.zeros(data.ndim)
    if upper is None:
        upper = data.shape
    lower = fix_limits(data, lower)
    upper = fix_limits(data, upper)
    m_slab = []
    for i in range(data.ndim):
        m_slab.append(slice(lower[i], upper[i]))
    return m_slab


def matching_slabs(data, flux, lower, upper):
    """
    Obtain the matching subcube inside the lower and upper points.

    Paramters
    ---------
    data : numpy.ndarray
        First data cube

    flux : numpy.ndarray
        Second data cubse

    lower : tuple
        Lower coordinates for the subcube.

    upper : tuple
        Upper coordinates for the subcube.


    Returns
    -------
    The subcube inside the lower and upper points that matches both data cube dimensions.
    """

    data_slab = slab(data, lower, upper)
    flow = np.zeros(flux.ndim)
    fup = np.array(flux.shape)
    for i in range(data.ndim):
        if data_slab[i].start == 0:
            flow[i] = flux.shape[i] - data_slab[i].stop
        if data_slab[i].stop == data.shape[i]:
            fup[i] = data_slab[i].stop - data_slab[i].start
    flux_slab = slab(flux, flow, fup)
    return data_slab, flux_slab

def index_mesh(data, lower=None, upper=None):
    """ Create an meshgrid from indices """
    sl = slab(data, lower, upper)
    dim = data.ndim
    slices = []
    for i in range(dim):
        slices.append(slice(sl[i].start, sl[i].stop))
    retval = np.mgrid[slices]
    return retval

def index_features(data, lower=None, upper=None):
    """ Creates an array with indices in features format """
    msh = index_mesh(data, lower, upper)
    dim = data.ndim
    ii = np.empty((dim, int(msh.size / dim)))
    for i in range(dim):
        ii[dim - i - 1] = msh[i].ravel()
    return ii