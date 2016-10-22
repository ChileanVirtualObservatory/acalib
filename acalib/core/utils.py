import numpy as np

def fix_mask(data, mask):
    ismasked = isinstance(data, np.ma.MaskedArray)
    if ismasked and mask is None:
        return data
    else:
        return np.ma.MaskedArray(data, mask)

def fix_limits(data, vect):
    """ Fix vect index to be inside data """
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
    """ Obtain the n-dimensional slab from lower to upper (i.e. slab is a vector of slices)"""
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
    """ Obtain the matching data and flux slabs from lower to upper while fixing the limits"""
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

