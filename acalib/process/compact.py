"""
.. module:: compact
   :platform: Unix
   :synopsis: compact representations of images

.. moduleauthor:: Mauricio Araya <maray@inf.utfsm.cl>

"""

import numpy as np
from collenctions import namedtuple

GMM=namedtuple('GMM',('amps','centers','pmat'))

def threshold_compact(data,threshold,cutoff):
    ff=np.where(data>threshold)
    amps=data[ff].filled()
    centers=np.transpose(ff).astype(float)
    I=np.identity(data.ndim)
    prec=[]
    for a in amps:
        val=2*np.log(a/cutoff)
        prec.append(val*I)
    prec=np.array(Prec)
    res=GMM(amps,centers,prec)
    return res

def isogaussian_compact(data):
    
    pass


