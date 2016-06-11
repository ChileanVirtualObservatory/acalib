"""
.. automodule:: gmr
   :platform: Unix
   :synopsis: Functions for ussing Gaussian Mixture Representations of data

.. moduleauthor:: Mauricio Araya <maray@inf.utfsm.cl>

"""

import numpy as np
from astropy.table import Table
from astropy import log
#from acalib.process.utils import *
from utils import *
import matplotlib.pyplot as plt


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

def _update_energies(energy,residual,mould,nlevel,delta,lower=None,upper=None):
    #def _update_energies(self,lb,ub):
    """Update the energies, only from the lower to the upper points. 
    """
    #TODO: I do now know if slab actually states that we are making a copy...
    if lower is None:
       lower=np.zeros(residual.ndim)
    lower=fix_limits(residual,lower)
    residual_slab=slab(residual,lower,upper)
    residual_view=residual[residual_slab]
    #resi_view=residual[energy_slab]
    feat=np.array(np.where(residual_view>nlevel))
    feat=feat.T
    #print mould
    #print lower,upper
    for idx in feat:
        base=idx + lower
        lb=base-delta
        ub=base+delta + 1
        #print base,lb,ub
        #print mould.shape
        residual_slab,mould_slab=matching_slabs(residual,mould,lb,ub)
        val=np.min(residual[residual_slab]/mould[mould_slab])
        #if val > 100:
        #    print("Base = "+str(base)+" Value = "+str(val))
        #    print residual[residual_slab]
        energy[tuple(base)]=val
        #energy[tuple(base)]=1000

def gmr_from_mould(data,threshold,nlevel,P,upper=None,lower=None,max_iter=None,full_output=False):
    debug=True
    inten=[]
    center=[]
    struct=[]
    niter=0
    message="Not finished"

    #Restrict data to what the corresponding slab
    data=data[slab(data,upper,lower)]
    if max_iter==None:
        d=data.ndim
        #GMRs are useless if we have more data than the original data
        max_iter=data.size/(1 + d + d*d)
    residual=data.copy()
    # Argmax
    max_idx=np.unravel_index(residual.argmax(),residual.shape)
    max_val=residual[max_idx]
    # Compute delta
    Sigma=np.linalg.inv(P)
    delta=np.sqrt(2*np.log(max_val/nlevel)*Sigma.diagonal())
    # Compute mould TODO
    mould=create_mould(P,delta)
    mould=mould/mould.max()
    #discretize delta
    delta=np.floor(np.array(mould.shape)/2)
    # Create energy matrix
    energy=np.zeros(data.shape)
    #np.empty_like(data)
    #mask=np.isnan(energy)
    #energy.data[np.logical_not(mask)]=max_val
    # check for diagonality
    #dstruct=np.all(a == np.diag(np.diag(P)))
    # compute 
    if debug:
       log.info("Delta = "+str(delta))
       log.info("Max Iter = "+str(max_iter))
       plt.imshow(mould.sum(axis=(0)))
       plt.show()
    _update_energies(energy,residual,mould,nlevel,delta)
    
    while True:
       niter+=1
       if (niter > max_iter):
           message="Maximum iterations reached = "+str(niter)
           break
       max_idx=np.array(np.unravel_index(energy.argmax(),energy.shape))
       max_val=energy[tuple(max_idx)]
       rem=max_val - nlevel
       
       if rem <= 0.0:
           message="No more signal to extract at iteration "+str(niter)
       if (rem < threshold):
           message="Signal reached the desired threshold of "+str(threshold)+" at iteration "+str(niter)
       inten.append(rem)
       center.append(max_idx)
       struct.append(P)
       lb=max_idx - delta
       ub=max_idx + delta +1 
       add_flux(residual,-rem*mould,lb,ub)
       _update_energies(energy,residual,mould,nlevel,delta,lower=lb-delta,upper=ub+delta)
       if debug:
           log.info("Iter "+str(niter))
           log.info("Maximum energy E = "+str(max_val)+" at "+str(max_idx))
           log.info("Remove E = "+str(rem)+" SNR = "+str(rem/nlevel))# + " GAP = "+ str(y/rms - 1.0 - snrlimit))
           print "max_energy = ",np.max(energy)
           print "min_energy = ",np.min(energy)
           print "max_resid = ",np.max(residual)
           print "min_resid = ",np.min(residual)
           plt.imshow(residual.sum(axis=(0)))
           plt.colorbar()
           plt.show()
           plt.imshow(energy.sum(axis=(0)))
           plt.colorbar()
           plt.show()
    inten=np.array(inten)
    center=np.array(center)
    struct=np.array(struct)
    res=Table([inten,center,struct],names=('intensity','center','struct'))
    if full_output:
       return res,message
    return res
  
def gmr_from_iterfit(data):
    pass

def gmr_from_heuristic(data):
    pass


if __name__ == '__main__':
    # SandBox space for testing
    a=np.random.random((200,200,200))
    P=np.array([[2,0,0],[0,1,0],[0,0,1]])
    freak=np.array([[2,0.5,0.9],[0.5,1,0.3],[0.9,0.3,1]])
    delta=np.array([50,50,50])
    p1=np.array([45,97,143])
    peak=create_mould(0.01*freak,delta)
    add_flux(a,10*peak,p1-delta,p1+delta+1)
    p2=np.array([143,45,97])
    peak=create_mould(0.05*freak,delta)
    add_flux(a,30*peak,p2-delta,p2+delta+1)
    plt.imshow(np.sum(a,axis=(0)))
    plt.colorbar()
    plt.show()
    (tab,message)=gmr_from_mould(a,3.0,1.0,0.5*P,full_output=True)
    print(message)
    print(tab)
