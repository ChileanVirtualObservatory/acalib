from astropy import log
from astropy.table import Table
import numpy as np
from . import matching_slabs, fix_limits, slab, snr_estimation, create_mould, add
from .function_models import _eighth_mould


def _update_min_energy(energy,mat,ub,lb,delta):
    """Updates the minimum energies of energy from mat defaced by delta. 
       ub and lb bounds are provided to shrink the mat matrix when required (out of bounds, or partial update)
    """
    # Numpyfy everithing
    ub=np.array(ub)
    lb=np.array(lb)
    delta=np.array(delta,dtype=int)
    # Create energy (e) and mat (m) indices 
    eub=ub + delta
    elb=lb + delta
    #print eub,elb,energy.shape
    eslab,mslab=matching_slabs(energy,mat,elb,eub)
    #mub=ub-lb
    #mub=mat.shape
    #mlb=np.array([0,0,0])
    #umask=eub > bord
    #lmask=elb < 0
    #mub[umask]-=eub[umask] - bord[umask]
    #mlb[lmask]-=elb[lmask]
    #eub[umask]=bord[umask]
    #elb[lmask]=0
    # Obtain a reduced view of the matrices
    eview=energy[eslab]
    mview=mat[mslab]
    #energy[elb[0]:eub[0],elb[1]:eub[1],elb[2]:eub[2]]
    #mview=mat[mlb[0]:mub[0],mlb[1]:mub[1],mlb[2]:mub[2]]
    #mview=mview.copy()
    # Select those that are lower in mat than in energy
    #print elb, eub, mlb, mub
    #print eview.shape,mview.shape
    # Update them in the energy matrix.
    #try:
    cmat=mview < eview
    #except ValueError:
    #   print energy.shape
    #   print elb,eub
    #   print eview.shape
    #   print mview.shape
    #   print mat.shape
    #   print eslab,mslab
    eview[cmat]=mview[cmat]


def _update_energies_sym(residual,energy,ev,ef,lb,ub):
      """Update the energies, only from the lb to the ub points. 
      """
      #TODO: I do now know if get_slice actually states that we are making a copy...
      lb=fix_limits(residual, lb)
      ub=fix_limits(residual, ub)
      mcb=residual[slab(residual, lb, ub)]
      #Obtain the reference of the eighth of the bubble.
      # Iterates for every point in the eighth of the bubble
      for i in range(0,ev.size):
         mat=mcb/ev[i]
         d=ef[i]
         dset=np.array([[1,1,1],[1,1,-1],[1,-1,1],[1,-1,-1],[-1,1,1],[-1,1,-1],[-1,-1,1],[-1,-1,-1]])*d
         dset=np.vstack({tuple(row) for row in dset})
         for delta in dset:
             _update_min_energy(energy,mat,ub,lb,delta)

def _precision_from_delta(delta,clev):
    delta=np.array(delta)
    sq_delta=1./(delta*delta)
    P=np.diag(sq_delta)
    return(-2*np.log(clev)*P)

def scatpix_detect(data,threshold=None,noise=None,upper=None,lower=None,full_output=False):
    """ Obtain an homogeneous representation using the scattered pixels over a threshold.

    This function generates an homogeneous representation by using only those pixels above the threshold. 
    Each pixel generates several identical values depending on the intensity of each pixel (i.e., floor(intensity/noise)). 

    :param data: n-dimensional array containing the data to be processed.
    :type data: ndarray
    :param threshold: the theshold to consider a pixel relevant to be included in the representation.
    :type threshold: float
    :param noise : noise level to be subtracted from the intensities 
    :type noise: float
    :returns: An astropy table with all the ::math:`\\mu` and the metadata.
    :rtype: Table

    """
# COMMENT Removed
#This can also be understood as very small Gaussians ::math:`G(x) = a \\exp(-0.5 (\mu - x)^\\top P (\mu - x))` 
#    where ::math:`a=\sigma` correspond to the noise parameter,  the center ::math:`\\mu` is the pixel position (pos) and ::math:`P` is a diagonal matrix of
#    the form ::math:`-2\\log(\delta) \\cdot I` (::math:`\delta` is a small number

    #Restrict data to what the corresponding slab
    data=data[slab(data,upper,lower)]

    ff=np.where(data>threshold)
    if isinstance(data,np.ma.MaskedArray):
        inten=data[ff].filled()
    else: 
        inten=data[ff]
    if full_output:
       residual=np.nan_to_num(data)
       synthetic=np.zeros(data.shape)
    ntimes=(inten/noise).astype(int)
    center=np.transpose(ff).astype(float)
    positions=[]
    for cen,tim in zip(center,ntimes):
        #print cen,tim,noise
        mylst=[cen.astype(int)]*tim
        positions.extend(mylst)
        if full_output:
            residual[tuple(cen.astype(int))]-=tim*noise
            synthetic[tuple(cen.astype(int))]=tim*noise
    positions=np.array(positions)
    rep=Table([positions],names=['center'])
    if full_output:
        return rep,synthetic,residual
    return rep

def bubble_detect(data,meta=None,noise=None,threshold=None,delta=None,gamma=0.1,full_output=False,verbose=False):
    if delta is None:
        if meta is None:
            delta=[1,1,1]
        else:
            spa=np.ceil((np.abs(meta['BMIN']/meta['CDELT1']) - 1)/2.0)
            delta=[1,spa,spa]
    if noise is None:
        noise=rms(data)
    if threshold is None:
        threshold=snr_estimation(data,mask=mask,noise=noise)*noise
    if verbose:
        print(threshold,noise,delta)
    P=_precision_from_delta(delta,gamma)
    mould=create_mould(P,delta)
    #equant=mould.sum()*noise
    residual=np.nan_to_num(data)
    energy=residual.copy()
    if full_output:
        synthetic=np.zeros(residual.shape)
        elist=[]
    (ev,ef)= _eighth_mould(P, delta)
    _update_energies_sym(residual,energy,ev,ef,lb=(0,0,0),ub=residual.shape)
    positions=[]
    niter=0
    delta=np.array(delta)
    while True:
        niter+=1
        idx      = np.unravel_index(energy.argmax(), energy.shape)
        max_ener = energy[idx]
        if verbose and niter%1000==0:
            log.info("Iteration: "+str(niter))
            log.info("Maximum energy E = "+str(max_ener)+" SNR = "+str(max_ener/noise))   #at "+str(idx))
            #fig = plt.figure(figsize=(20,5))
            #ax = fig.add_subplot(121)
            #ax.imshow(integrate(residual).data, origin='lower')
            #ax.set_title('Residual')
            #ax = fig.add_subplot(122)
            #ax.imshow(integrate(synthetic).data, origin='lower')
            #ax.set_title('Denoised')
            #plt.show()
       
        if max_ener < noise:
            if verbose:
                log.info("Criterion Met: Energy < Noise Level ")
            break
        if (max_ener < threshold):
            if verbose:
                log.info("Criterion Met: SNR="+str(max_ener/noise)+"<"+str(threshold/noise))
            break
        ub=idx + delta + 1
        lb=idx - delta
        add(residual,-noise*mould,lb,ub)
        if full_output:
            add(synthetic,noise*mould,lb,ub)
            elist.append(max_ener)
        _update_energies_sym(residual,energy,ev,ef,lb,ub)
        positions.append(idx)
    positions=np.array(positions)
    rep=Table([positions],names=['center'])
    if full_output:
        return rep,synthetic,residual,energy,elist
    return rep

    # DO THIS FOR BUBBLE
def synthesize_bubbles(syn,pos,mould,nlevel,delta):
    for idx in pos:
        ub=idx + delta + 1
        lb=idx - delta
        add(syn,nlevel*mould,lb,ub)
    return syn

def precision_from_delta(delta,clev):
    delta=np.array(delta)
    sq_delta=1./(delta*delta)
    P=np.diag(sq_delta)
    return(-2*np.log(clev)*P)

