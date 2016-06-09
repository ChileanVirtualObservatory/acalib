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
    """ Obtain a pixel-based Gaussian Mixture Representation (GMR) from data.

    Args:
       data (ndarray): n-dimensional array containing the data to be processed.
       threshold (float): the theshold to consider a pixel relevant to be included in the representation.
       nlevel (float): noise level to be subtracted from the intensities and to compute the shape of each Gaussian.

    Returns:
       Table. An astropy table of m elements with the columns 
         * intensity (float)
         * center (ndarray) of size 
         * shape

    This function generate a GMR by using only those pixels above the threshold. Each pixel generates a very small
    Gaussian:

    $ G(x) = a \exp(-0.5 (\mu - x)^\top P (\mu - x)) $

    where $a$ is the instensity of the pixel minus nlevel, the center $\mu$ is the pixel position and $P$ is a diagonal matrix of
    the form:
    
    $ $ 

    """    
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
