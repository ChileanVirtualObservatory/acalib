
import numpy as np
from astropy.table import Table
from astropy import log
from astropy import units as u
from utils import *
import matplotlib.pyplot as plt
from convert import *

def _eighth_mould(P,delta):
    """This function creates a Gaussian mould with precision matrix P, using the already computed values of delta
    """
    n=len(delta)
    ax=[]
    elms=[]
    for i in range(n):
        lin=np.linspace(0,delta[i],delta[i]+1)
        elms.append(len(lin))
        ax.append(lin)
    grid=np.meshgrid(*ax,indexing='ij')
    feat=np.empty((n,np.product(elms)))
    for i in range(n):
        feat[i]=grid[i].ravel()
    mould=gaussian_function(np.zeros(n),P,feat,1)
    return mould,feat.T

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
    try:
       cmat=mview < eview
    except ValueError:
       print energy.shape
       print elb,eub
       print eview.shape
       print mview.shape
       print mat.shape
       print eslab,mslab
    eview[cmat]=mview[cmat]


def _update_energies_sym(residual,energy,ev,ef,lb,ub):
      """Update the energies, only from the lb to the ub points. 
      """
      #TODO: I do now know if get_slice actually states that we are making a copy...
      lb=fix_limits(residual,lb)
      ub=fix_limits(residual,ub)
      mcb=residual[slab(residual,lb,ub)]
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

@support_nddata
def scatpix_detect(data,wcs=None,meta=None,threshold=None,noise=None,upper=None,lower=None,full_output=False):
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

@support_nddata
def bubble_detect(data,wcs=None,meta=None,noise=None,threshold=None,delta=None,gamma=0.1,full_output=False,verbose=False):
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
        print threshold,noise,delta
    P=_precision_from_delta(delta,gamma)
    mould=create_mould(P,delta)
    #equant=mould.sum()*noise
    residual=np.nan_to_num(data)
    energy=residual.copy()
    if full_output:
        synthetic=np.zeros(residual.shape)
        elist=[]
    (ev,ef)=_eighth_mould(P,delta)
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
        add_flux(residual,-noise*mould,lb,ub)
        if full_output:
            add_flux(synthetic,noise*mould,lb,ub)
            elist.append(max_ener)
        _update_energies_sym(residual,energy,ev,ef,lb,ub)
        positions.append(idx)
    positions=np.array(positions)
    rep=Table([positions],names=['center'])
    if full_output:
        return rep,synthetic,residual,energy,elist
    return rep

    


#def gmr_from_mould(data,threshold,nlevel,P,upper=None,lower=None,max_iter=None,full_output=False):
#    debug=True
#    inten=[]
#    center=[]
#    struct=[]
#    niter=0
#    message="Not finished"
#
#    #Restrict data to what the corresponding slab
#    data=data[slab(data,upper,lower)]
#    if max_iter==None:
#        d=data.ndim
#        #GMRs are useless if we have more data than the original data
#        max_iter=data.size/(1 + d + d*d)
#    residual=data.copy()
#    # Argmax
#    max_idx=np.unravel_index(residual.argmax(),residual.shape)
#    max_val=residual[max_idx]
#    # Compute delta
#    Sigma=np.linalg.inv(P)
#    delta=np.round(np.sqrt(2*np.log(max_val/nlevel)*Sigma.diagonal()))
#    # Compute mould TODO
#    mould=create_mould(P,delta)
#    #plt.imshow(mould[0,:,:])
#    #plt.colorbar()
#    #plt.show()
#    #plt.imshow(mould[1,:,:])
#    #plt.colorbar()
#    #plt.show()
#    #plt.imshow(mould[2,:,:])
#    #plt.colorbar()
#    #plt.show()
#    #plt.imshow(mould[3,:,:])
#    #plt.colorbar()
#    #plt.show()
#    #plt.imshow(mould[4,:,:])
#    #plt.colorbar()
#    #plt.show()
#    #mould=mould/mould.max()
#    #discretize delta
#    delta=np.floor(np.array(mould.shape)/2)
#    # Create energy matrix
#    energy=np.zeros(data.shape)
#    #np.empty_like(data)
#    #mask=np.isnan(energy)
#    #energy.data[np.logical_not(mask)]=max_val
#    # check for diagonality
#    #dstruct=np.all(a == np.diag(np.diag(P)))
#    # compute 
#    if debug:
#       log.info("Delta = "+str(delta))
#       log.info("Max Iter = "+str(max_iter))
#       #plt.imshow(mould.sum(axis=(0)))
#       #plt.show()
#    _update_energies(energy,residual,mould,nlevel,delta)
#    
#    while True:
#       niter+=1
#       if (niter > max_iter):
#           message="Maximum iterations reached = "+str(niter)
#           break
#       max_idx=np.array(np.unravel_index(energy.argmax(),energy.shape))
#       max_val=energy[tuple(max_idx)]
#       #    tlb=base-delta
#       #    tub=base+delta + 1
#       #    plt.imshow(residual[slab(residual,tlb,tub)].sum(axis=0))
#       #    print energy[tuple(base)]
#       #    print residual[slab(residual,tlb,tub)]
#       #    print mould
#       #    print (residual[slab(residual,tlb,tub)]/mould)
#       #    plt.show()
#       rem=max_val - nlevel
#       
#       if rem <= 0.0:
#           message="No more signal to extract at iteration "+str(niter)
#           break
#       if (rem < threshold):
#           message="Signal reached the desired threshold of "+str(threshold)+" at iteration "+str(niter)
#           break
#       inten.append(rem)
#       center.append(max_idx)
#       struct.append(P)
#       lb=max_idx - delta
#       ub=max_idx + delta +1 
#       add_flux(residual,-rem*mould,lb,ub)
#       _update_energies(energy,residual,mould,nlevel,delta,lower=lb-delta,upper=ub+delta)
#       if debug:
#           log.info("Iter "+str(niter))
#           log.info("Maximum energy E = "+str(max_val)+" at "+str(max_idx))
#           log.info("Remove E = "+str(rem)+" SNR = "+str(rem/nlevel))# + " GAP = "+ str(y/rms - 1.0 - snrlimit))
#           print "max_energy = ",np.max(energy)
#           print "min_energy = ",np.min(energy)
#           print "max_resid = ",np.max(residual)
#           print "min_resid = ",np.min(residual)
#           base=np.unravel_index(residual.argmax(),residual.shape)
#           print "e(max_resid)=",energy[tuple(base)]
#           #plt.imshow(residual.sum(axis=(0)))
#           #plt.colorbar()
#           #plt.show()
#           #plt.imshow(energy.sum(axis=(0)))
#           #plt.colorbar()
#           #plt.show()
#    inten=np.array(inten)
#    center=np.array(center)
#    struct=np.array(struct)
#    res=Table([inten,center,struct],names=('intensity','center','struct'))
#    if full_output:
#       return res,message,residual,energy
#    return res
# 
# 
#
#if __name__ == '__main__':
#    # SandBox space for testing
#    a=3*np.random.random((50,200,200))
#    P=np.array([[2,0,0],[0,1,0],[0,0,1]])
#    freak=np.array([[2,0.5,0.9],[0.5,1,0.3],[0.9,0.3,1]])
#    delta=np.array([10,50,50])
#    p1=np.array([15,97,143])
#    peak=create_mould(0.01*freak,delta)
#    add_flux(a,5*peak,p1-delta,p1+delta+1)
#    p2=np.array([30,45,97])
#    peak=create_mould(0.02*freak,delta)
#    add_flux(a,10*peak,p2-delta,p2+delta+1)
#    rms=estimate_rms(a)
#    print("RMS = "+str(rms))
#    plt.imshow(np.sum(a,axis=(0)))
#    plt.colorbar()
#    plt.show()
#    (tab,message,residual,energy)=gmr_from_mould(a,0.2*rms,rms,P,full_output=True)
#    print(message)
#    print(tab)
#    plt.imshow(np.sum(residual,axis=(0)))
#    plt.colorbar()
#    plt.show()
#    plt.imshow(np.sum(energy,axis=(0)))
#    plt.colorbar()
#    plt.show()
#
#
