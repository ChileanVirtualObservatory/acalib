import numpy as np
import matplotlib.pyplot as plt
from astropy import log
import astropy.units as u
from astropy.nddata import *
import scipy.ndimage.interpolation as sni

def fix_mask(data,mask):
    ismasked=isinstance(data,np.ma.MaskedArray)
    if ismasked and mask is None: 
        return data
    else:
       return np.ma.MaskedArray(data,mask)     

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
#         print data_slab
#          print flux.shape
          flow[i] = flux.shape[i] - data_slab[i].stop
       if data_slab[i].stop == data.shape[i]:
#          print data_slab
#          print flux.shape
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

@support_nddata
def world_features(data,wcs,lower=None,upper=None):
    ii=to_features(data,lower,upper)
    f=wcs.wcs_pix2world(ii.T,0)
    f=f.T
    return f

