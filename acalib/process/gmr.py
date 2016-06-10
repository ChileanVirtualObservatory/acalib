"""
.. automodule:: gmr
   :platform: Unix
   :synopsis: Functions for ussing Gaussian Mixture Representations of data

.. moduleauthor:: Mauricio Araya <maray@inf.utfsm.cl>

"""

import numpy as np
from astropy.table import Table
from astropy import log
from acalib.process.utils import *

def gmr_from_pixels(data,threshold,nlevel,upper=None,lower=None):
    """ Obtain a pixel-based Gaussian Mixture Representation (GMR) from data.

    This function generate a GMR by using only those pixels above the threshold. Each pixel generates a very small
    Gaussian ::math:`G(x) = a \\exp(-0.5 (\mu - x)^\\top P (\mu - x))` 
    where ::math:`a` is the instensity (data[pos]) of the pixel ::math:`n_{level}` (nlevel), the center ::math:`\\mu` is the pixel position (pos) and ::math:`P` is a diagonal matrix of
    the form ::math:`2\\log(a/n_{level}) \\cdot I`

    :param data: n-dimensional array containing the data to be processed.
    :type data: ndarray
    :param threshold: the theshold to consider a pixel relevant to be included in the representation.
    :type threshold: float
    :param nlevel: noise level to be subtracted from the intensities and to compute the structure of each Gaussian.
    :type nlevel: float
    :param upper: 
    :type upper: ndarray
    :returns: An astropy table, where each row has the parameters of a single Gaussian ::math:`a`, ::math:`\\mu` and ::math:`P` (intensity, center and structure).
    :rtype: Table

    :Example:
    
    >>> a=np.random.random((2,2,2))
    >>> tab=gmr_from_pixels(a,0.5,0.001)
    intensity      center [3]        structure [3,3]          
    -------------- ---------- ------------------------------
    0.537338435568 0.0 .. 0.0 12.5732562597 .. 12.5732562597
    0.685661594975 0.0 .. 0.0 13.0607684085 .. 13.0607684085
    0.939673095857 0.0 .. 1.0 13.6910640888 .. 13.6910640888
    0.554681589695 1.0 .. 1.0 12.6367884737 .. 12.6367884737
    0.522312859713 1.0 .. 0.0 12.5165335129 .. 12.5165335129

    """
    #Restrict data to what the corresponding slab
    data=data[slab(data,upper,lower)]

    ff=np.where(data>threshold)
    if isinstance(data,np.ma.MaskedArray):
        inten=data[ff].filled()
    else: 
        inten=data[ff]
    inten-=nlevel
    center=np.transpose(ff).astype(float)
    I=np.identity(data.ndim)
    struct=[]
    for i in inten:
        val=2*np.log(i/nlevel)
        struct.append(val*I)
    struct=np.array(struct)
    res=Table([inten,center,struct],names=('intensity','center','structure'))
    return res

#_compute_bubble   
#_update_energies(energy,bub)

def gmr_from_mould(data,threshold,nlevel,P,verbose=False,upper=None,lower=None,max_iter=None,full_output=False):
    debug=False
    inten=[]
    center=[]
    struct=[]
    niter=0
    message="Not finished"

    #Restrict data to what the corresponding slab
    data=data[slab(data,upper,lower)]
    if max_iter=None
        max_iter

    if max_iter==None:
    residual=data.copy()
    energy=data.copy()
    # Argmax
    max_idx=np.unravel_index(residual.argmax(),residual.shape)
    max_val=residual[max_idx]
    # Compute delta
    delta=np.sqrt(2*np.log(datamax/nlevel)*P.diagonal())))
    # Compute bubble TODO
    bub=_compute_bubble()
    # Create energy matrix
    mask=np.isnan(energy)
    energy.data[np.logical_not(mask)]=max_val
    _update_energies(energy,bub)
    
    while True:
       niter+=1
       if (niter > max_iter):
           message="Maximum iterations reached = "+str(maxbub))
           break
       max_idx=np.array(np.unravel_index(residual.argmax(),residual.shape))
       max_val=residual[max_idx]
       xmax=np.array(xmax)
       rem=max_val - nlevel
       if rem <= 0.0:
           if verbose:
               message="No more signal to extract at iteration "+str(niter)
           break
       if (rem < threshold):
           if verbose:
               message="Signal reached the desired threshold of "+str(threshold)+" at iteration "+str(niter)
           break
       inten.append(rem)
       center.append(xmax)
       struct.append(P)
       if debug:
           log.info("Iter "+str(niter))
           log.info("Maximum energy E = "+str(max_val)+" at "+str(xmax))
           log.info("Remove E = "+str(rem)+" SNR = "+str(rem/nlevel)# + " GAP = "+ str(y/rms - 1.0 - snrlimit))
       ub=xmax + delta + 1
       lb=xmax - delta
       add_flux(residual,-rem*bub,lb,ub)
       _update_energies(energy,lb,ub)
    inten=np.array(inten)
    center=np.array(center)
    struct=np.array(struct)
    res=Table([inten,center,struct],names=('intensity','center','struct'))
    return res
  
def gmr_from_iterfit(data):
    pass

def gmr_from_heuristic(data):
    pass


if __name__ == '__main__':
    # SandBox space for testing
    a=np.random.random((20,20,20))
    tab=gmr_from_mould(a,0.2,0.01,eye(3))
    print(tab)
