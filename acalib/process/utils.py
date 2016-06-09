"""
.. automodule:: utils
   :platform: Unix
   :synopsis: Helper functions for processing data

.. moduleauthor:: Mauricio Araya <maray@inf.utfsm.cl>

"""

import numpy as np

def fix_limits(data,vect):
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


def slab(data,lower,upper):
    if lower==None:
        lower=np.zeros(data.ndim)
    if upper==None:
        upper=self.data.shape
    lower=fix_limits(data,lower)
    upper=fix_limits(data,upper)
    m_slab=[]
    for i in range(data.ndim):
       m_slab.append(slice(lower[i],upper[i]))
    return m_slab

# TODO Throw exceptions...
def add_flux(data,flux,lower=None,upper=None):
    #if data.ndim!=flux.ndim:
    #    log.error("")
    data_slab=slab(data,lower,upper)
    flux_slab=slab(flux)
    for i in range(data.ndim):
       if data_slab[i].start == 0:
          flux_slab[i].start = flux.shape[i] - data_slab[i].stop
       if data_slab[i].stop == data.shape[i]:
          flux_slab[i].stop = data_slab[i].stop - data_slab[i].start
    self.data[slab]+=flux[flux_slab]

