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

   def _update_min_energy(self,mat,ub,lb,delta):
     """Updates the minimum energies of self.energy from mat defaced by delta. 
        ub and lb bounds are provided to shrink the mat matrix when required (out of bounds, or partial update)
     """
     bord=np.array(self.energy.shape())
     # TODO: Bad usage of AData, because we want a reference of the data to modify it
     ene=self.energy.data
     # Numpyfy everithing
     ub=np.array(ub)
     lb=np.array(lb)
     delta=np.array(delta,dtype=int)
     # Create energy (e) and mat (m) indices 
     eub=ub + delta
     elb=lb + delta
     mub=ub-lb
     mlb=np.array([0,0,0])
     umask=eub > bord
     lmask=elb < 0
     mub[umask]-=eub[umask] - bord[umask]
     mlb[lmask]-=elb[lmask]
     eub[umask]=bord[umask]
     elb[lmask]=0
     # Obtain a reduced view of the matrices
     eview=ene[elb[0]:eub[0],elb[1]:eub[1],elb[2]:eub[2]]
     mview=mat[mlb[0]:mub[0],mlb[1]:mub[1],mlb[2]:mub[2]]
     mview=mview.copy()
     # Select those that are lower in mat than in energy
     #print elb, eub, mlb, mub
     #print eview.shape,eview.shape
     cmat=mview < eview
     #print cmat
     #print eview[cmat]
     #print mview
     # Update them in the energy matrix.
     try:
        a=mview[cmat]
     except IndexError:
        print eview.shape
        print mview.shape
        print mat.shape
        print elb,eub,mlb,mub
     eview[cmat]=mview[cmat]


def _update_energies(energy,mould,lower=None,upper=None,dstruct=False):
      #def _update_energies(self,lb,ub):
      """Update the energies, only from the lb to the ub points. 
      """
      #TODO: I do now know if get_slice actually states that we are making a copy...
      lb=self.residual.fix_limits(lb)
      ub=self.residual.fix_limits(ub)
      mcb=self.residual.cut(lb,ub)
      #Obtain the reference of the eighth of the bubble.
      vv=self.eival
      ff=self.eifeat.T
      # Iterates for every point in the eighth of the bubble
      # this starts from one because we do not want to repeat the position (0,0,0)
      delta=np.array([1,1,1])*0
      mat=mcb/vv[0]
      self._update_min_energy(mat,ub,lb,delta)
      for i in range(1,vv.size):
         mat=mcb/vv[i]
         d=ff[i]
         # update in the eight directions
         delta=np.array([1,1,1])*d
         self._update_min_energy(mat,ub,lb,delta)
         delta=np.array([1,1,-1])*d
         self._update_min_energy(mat,ub,lb,delta)
         delta=np.array([1,-1,1])*d
         self._update_min_energy(mat,ub,lb,delta)
         delta=np.array([1,-1,-1])*d
         self._update_min_energy(mat,ub,lb,delta)
         delta=np.array([-1,1,1])*d
         self._update_min_energy(mat,ub,lb,delta)
         delta=np.array([-1,1,-1])*d
         self._update_min_energy(mat,ub,lb,delta)
         delta=np.array([-1,-1,1])*d
         self._update_min_energy(mat,ub,lb,delta)
         delta=np.array([-1,-1,-1])*d
         self._update_min_energy(mat,ub,lb,delta)

    

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
    # Compute mould TODO
    mould=create_mould(P,delta)
    # Create energy matrix
    mask=np.isnan(energy)
    energy.data[np.logical_not(mask)]=max_val
    # check for diagonality
    dstruct=np.all(a == np.diag(np.diag(P)))
    # compute 
    _update_energies(energy,mould,dstruct=dstruct)
    
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
       add_flux(residual,-rem*bub,lb,ub,dstruct=dstruct)
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
