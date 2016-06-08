"""
.. module:: gmr
   :platform: Unix
   :synopsis: Functions for ussing Gaussian Mixture Representations of data

.. moduleauthor:: Mauricio Araya <maray@inf.utfsm.cl>

"""

import numpy as np
#from collenctions import namedtuple
from astropy.table import Table


def gmr_from_pixels(data,threshold,nlevel):
    
    ff=np.where(data>threshold)
    if isinstance(data,np.ma.MaskedArray):
        inten=data[ff].filled()
    else: 
        inten=data[ff]
    intensity-=nlevel
    center=np.transpose(ff).astype(float)
    I=np.identity(data.ndim)
    shape=[]
    for i in inten:
        val=2*np.log(i/nlevel)
        shape.append(val*I)
    shape=np.array(shape)
    res=Table([intensity,center,shape],names=('intensity','center','shape'))
    return res

def gmr_from_bubbles(data,threshold,nlevel):
    pass

def gmr_from_iterfit(data):
    pass

def gmr_from_heuristic(data):


if __name__ == '__main__':
    # SandBox space for testing
    a=np.random.random((10,20,300))
    tab=threshold_compact(a,0.1,0.001)
