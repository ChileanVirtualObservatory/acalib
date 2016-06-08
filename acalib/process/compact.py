"""
.. module:: compact
   :platform: Unix
   :synopsis: compact representations of images

.. moduleauthor:: Mauricio Araya <maray@inf.utfsm.cl>

"""

import numpy as np
#from collenctions import namedtuple
from astropy.table import Table

def threshold_compact(data,threshold,cutoff):
    ff=np.where(data>threshold)
    if isinstance(data,np.ma.MaskedArray):
        amps=data[ff].filled()
    else: 
        amps=data[ff]
    centers=np.transpose(ff).astype(float)
    I=np.identity(data.ndim)
    prec=[]
    for a in amps:
        val=2*np.log(a/cutoff)
        prec.append(val*I)
    prec=np.array(prec)
    res=Table([amps,centers,prec],names=('amps','centers','prec'))
    return res

def isogaussian_compact(data):
    pass

if __name__ == '__main__':
    # SandBox space for testing
    a=np.random.random((10,20,300))
    tab=threshold_compact(a,0.1,0.001)
