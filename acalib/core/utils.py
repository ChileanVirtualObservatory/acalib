import numpy as np

def fix_mask(data, mask):
    """

    Parameters
    ----------
    data: numpy.ndarray or numpy.ma.MaskedArray
        Astronomical data cube.
    mask: numpy.ndarray
        Boolean that will be applied.
    
    Returns
    -------
    result: numpy.ma.MaskedArray
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
    data: numpy.ndarray or numpy.ma.MaskedArray
        Astronomical data cube.
    vect: tuple, list or numpy.ndarray
        Array with the indexes to be fixed.

    Returns
    -------
    result: numpy.ndarray
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
    data: numpy.ndarray
        Atronomical data cube.
    lower: 3-tuple (default=None)
    upper: 3-tuple (default=None)


    Returns
    -------
    result: list
       Sub-cube.
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
    Obtain the matching data and flux slabs from lower to upper while fixing the limits
    
    Paramters
    ---------
    data: numpy.ndarray
    flux: numpy.ndarray
    lower: tuple
    upper: tuple


    Returns
    -------
    result: tuple
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

