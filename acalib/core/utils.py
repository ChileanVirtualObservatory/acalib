
import numpy as np
import matplotlib.pyplot as plt
from indices import *
import astropy.units as u

def add_flux(data,flux,lower=None,upper=None):
    """ Adds flux to data. 

    Lower and upper are bounds for data. This operation is border-safe. 
    """
    #if data.ndim!=flux.ndim:
    #    log.error("")
    data_slab,flux_slab=matching_slabs(data,flux,lower,upper)
    data[data_slab]+=flux[flux_slab]

def gaussian_function(mu,P,feat,peak):
    """ Generates an n-dimensional Gaussian using the feature matrix feat,
    centered at mu, with precision matrix P and with intensity peak.
    """
    #print feat
    cent_feat=np.empty_like(feat)
    for i in range(len(mu)):
       cent_feat[i]=feat[i] - mu[i]
    qform=(P.dot(cent_feat))*cent_feat
    quad=qform.sum(axis=0)
    res=np.exp(-quad/2.0)
    res=peak*(res/res.max())
    return res

# TODO: extend to n-dimensions (only works for 3)
def axes_ranges(data,wcs,lower=None,upper=None):
    """ Get axes extent (transforms freq to velocity!) """
    if lower==None:
        lower=[0,0,0]
    if upper==None:
        upper=data.shape
    lower=lower[::-1]
    lwcs=wcs.wcs_pix2world([lower], 0)
    lwcs=lwcs[0]
    upper=upper[::-1]
    uwcs=wcs.wcs_pix2world([upper], 0)
    uwcs=uwcs[0]
    lfreq=lwcs[2]*u.Hz
    ufreq=uwcs[2]*u.Hz
    rfreq=wcs.wcs.restfrq*u.Hz
    eq= u.doppler_radio(rfreq)
    lvel=lfreq.to(u.km/u.s, equivalencies=eq)
    uvel=ufreq.to(u.km/u.s, equivalencies=eq)
    ranges=[lvel.value,uvel.value,lwcs[1],uwcs[1],lwcs[0],uwcs[0]]
    return ranges

#TODO: try to merge with axes_ranges!
def axis_range(data,wcs,axis):
    lower=wcs.wcs_pix2world([[0,0,0]], 0) - wcs.wcs.cdelt/2.0
    shape=data.shape
    shape=[shape[::-1]]
    upper=wcs.wcs_pix2world(shape, 1) + wcs.wcs.cdelt/2.0
    return (lower[0][axis],upper[0][axis])  


def create_mould(P,delta):
    """This function creates a Gaussian mould with precision matrix P, using the already computed values of delta
    """
    n=len(delta)
    ax=[]
    elms=[]
    for i in range(n):
        lin=np.linspace(-delta[i]-0.5,delta[i]+0.5,delta[i]*2+1)
        elms.append(len(lin))
        ax.append(lin)
    grid=np.meshgrid(*ax,indexing='ij')
    feat=np.empty((n,np.product(elms)))
    for i in range(n):
        feat[i]=grid[i].ravel()
    mould=gaussian_function(np.zeros(n),P,feat,1)
    mould=mould.reshape(*elms)
    return(mould)

def estimate_rms(data):
    """A simple estimation of the RMS
    """
    mm=data * data
    if isinstance(mm,np.ma.MaskedArray):
        rms=np.sqrt(mm.sum()*1.0/mm.count())
    else:
        rms=np.sqrt(mm.sum()*1.0/mm.size)
    return rms

def gaussflux_from_world_window(data,wcs,mu,P,peak,cutoff):
   Sigma=np.linalg.inv(P)
   window=np.sqrt(2*np.log(peak/cutoff)*np.diag(Sigma))
   lower,upper=world_window_to_index(data,wcs,mu,window)
   if np.any(np.array(upper-lower)<=0):
       return None,lower,upper
   feat=world_features(data,wcs,lower,upper)
   res=gaussian_function(mu,P,feat,peak)
   res=res.reshape(upper[0]-lower[0],upper[1]-lower[1],upper[2]-lower[2])
   return res,lower,upper

def world_features(data,wcs,lower=None,upper=None):
    ii=to_features(data,lower,upper)
    f=wcs.wcs_pix2world(ii.T,0)
    f=f.T
    return f


#if __name__ == '__main__':
#    # Slab and AddFlux test
#    a=np.random.random((20,20,20))
#    sl=slab(a,(-5,4,5),(15,25,10))
#    print(sl)
#    b=100*np.random.random((10,10,10))
#    add_flux(a,b,(15,-5,7),(25,5,17))
#    c=np.where(a>1.0)
#    print(str(c[0].size)+" should be near 250")
#    # Mould test
#    P=np.array([[0.05,0.01,0],[0.01,0.07,0.03],[0,0.03,0.09]])
#    delta=[10,15,20]
#    mould=create_mould(P,delta)
#    plt.imshow(mould.sum(axis=(0)))
#    plt.show()




