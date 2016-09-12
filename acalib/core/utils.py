import numpy as np
import matplotlib.pyplot as plt
from indices import *
from astropy import log
import astropy.units as u
from astropy.nddata import *
import scipy.ndimage.interpolation as sni

def fix_mask(data,mask):
    ismasked=isinstance(data,np.ma.MaskedArray)
    if ismasked and mask is None: 
        return data
    else:
       return np.ma.MaskedArray(data,mask)     

@support_nddata
def rotate(data,angle):
    return sni.rotate(data,angle)

@support_nddata
def standarize(data,wcs=None,unit=None,mask=None,meta=None):
    if mask is not None:
        data=fix_mask(data,mask)
    y_min=data.min()
    res=data - y_min
    y_fact=res.sum()
    res=res/y_fact
    nres=NDData(res, uncertainty=None, mask=mask,wcs=wcs, meta=meta, unit=unit)
    return (nres,y_min,y_fact)

# TODO need to be nddatafied
def unstandarize(data, y_min,y_fact):
    return data*y_fact + y_min

@support_nddata
def cut(data,wcs=None,mask=None,unit=None,lower=None,upper=None):
    mslab=slab(data,lower,upper)
    scube=data[mslab]
    newwcs=wcs.slice(mslab,numpy_order=True)
    return NDData(scube,wcs=newwcs,unit=unit)

# TODO: generalize this function... is not very generic :S
@support_nddata
def spectra(data,wcs=None,mask=None,unit=None,position=None,aperture=None):
    if position is None:
        # Get celestial center
        position=wcs.celestial.wcs.crval*u.deg
    if aperture is None:
        # Get 1 pixel aperture
        aperture=np.abs(wcs.celestial.wcs.cdelt[0])*u.deg 
    if position.unit == u.pix and aperture.unit == u.pix:
        # TODO:  Here is the nasty part
        lb=np.array([0,            position[1].value - aperture.value, position[0].value - aperture.value])
        ub=np.array([data.shape[2],position[1].value + aperture.value, position[0].value + aperture.value])
    else:
        log.error("Not Implemented Yet!")
    specview=data[slab(data,lb,ub)]
    return specview.sum(axis=(1,2))

def moment(data,order,wcs=None,mask=None,unit=None,restfrq=None):
    if wcs is None:
        log.error("A world coordinate system (WCS) is needed")
        return None
    data=fix_mask(data,mask)
    dim=wcs.wcs.spec
    rdim=data.ndim - 1 - dim
    v=get_velocities(data,wcs,np.arange(data.shape[rdim]),restfrq)
    v=v.value
    #delta=np.mean(np.abs(v[:v.size-1] - v[1:v.size]))
    #newdata=data.sum(axis=rdim)*delta
    m0=data.sum(axis=rdim)
    if order==0:
        mywcs=wcs.dropaxis(dim)
        return NDData(m0.data, uncertainty=None, mask=m0.mask,wcs=mywcs, meta=None, unit=unit)
    #mu,alpha=np.average(data,axis=rdim,weights=v,returned=True)
    mu,alpha=np.ma.average(data,axis=rdim,weights=v,returned=True)
    m1=alpha*mu/m0
    if order==1:
        mywcs=wcs.dropaxis(dim)
        return NDData(m1.data, uncertainty=None, mask=m1.mask,wcs=mywcs, meta=None, unit=u.km/u.s)
    v2=v*v
    var,beta=np.ma.average(data,axis=rdim,weights=v2,returned=True)
    #var,beta=data.average(axis=rdim,weights=v2,returned=True)
    m2=np.sqrt(beta*var/m0 - m1*m1)
    if order==2:
        mywcs=wcs.dropaxis(dim)
        return NDData(m2.data, uncertainty=None, mask=m2.mask,wcs=mywcs, meta=None, unit=u.km*u.km/u.s/u.s)
    log.error("Order not supported")
    return None
        
# Should return a NDData
@support_nddata
def moment0(data,wcs=None,mask=None,unit=None,restfrq=None):
    return moment(data,0,wcs,mask,unit,restfrq)

@support_nddata
def moment1(data,wcs=None,mask=None,unit=None,restfrq=None):
    return moment(data,1,wcs,mask,unit,restfrq)

@support_nddata
def moment2(data,wcs=None,mask=None,unit=None,restfrq=None):
    return moment(data,2,wcs,mask,unit,restfrq)


@support_nddata
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

@support_nddata
def denoise(data,wcs=None,mask=None,unit=None,threshold=0.0):
      elms=data>threshold
      newdata=np.zeros(data.shape)
      newdata[elms]=data[elms]
      return NDData(newdata, uncertainty=None, mask=mask,wcs=wcs, meta=None, unit=unit)

@support_nddata
def get_velocities(data,wcs=None,fqi=None,restfrq=None):
    if wcs is None:
        log.error("A world coordinate system (WCS) is needed")
        return None
    if fqi is None:
        return None
    if restfrq is None:
        restfrq=wcs.wcs.restfrq*u.Hz
    dim=wcs.wcs.spec
    idx=np.zeros((fqi.size,data.ndim))
    idx[:,dim]=fqi
    vals=wcs.all_pix2world(idx,0)
    eq=u.doppler_radio(restfrq)
    vec=vals[:,dim]*u.Hz
    return vec.to(u.km/u.s, equivalencies=eq)

# TODO: extend to n-dimensions (only works for 3)
@support_nddata
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

#TODO: try to merge with axes_ranges and get_velocities!
@support_nddata
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
        lin=np.linspace(-delta[i],delta[i],delta[i]*2+1)
        elms.append(len(lin))
        ax.append(lin)
    grid=np.meshgrid(*ax,indexing='ij')
    feat=np.empty((n,np.product(elms)))
    for i in range(n):
        feat[i]=grid[i].ravel()
    mould=gaussian_function(np.zeros(n),P,feat,1)
    mould=mould.reshape(*elms)
    return mould


@support_nddata
def rms(data,mask=None):
    """A simple estimation of the noise level by computing the RMS. If mask != None, then 
       we use that mask.
    """
    if mask is not None:
        data=fix_mask(data,mask)
    mm=data*data
    #if mask is not None and not ismasked:
    rms=np.sqrt(mm.sum()*1.0/mm.size)
    return rms

@support_nddata
def snr_estimation(data,mask=None,noise=None,points=1000,full_output=False):
    if noise is None:
       noise=rms(data,mask)
    x=[]
    y=[]
    sdata=data[data>noise]
    for i in range(1,int(points)):
        val=1.0 + 2.0*i/points
        sdata=sdata[sdata>val*noise]
        if sdata.size < 2:
            break
        yval=sdata.mean()/noise
        x.append(val)
        y.append(yval)
    y=np.array(y)
    v=y[1:]-y[0:-1]
    p=v.argmax() + 1
    snrlimit=x[p]
    if full_output==True:
       return snrlimit,noise,x,y,v,p 
    return snrlimit

@support_nddata
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

@support_nddata
def world_features(data,wcs,lower=None,upper=None):
    ii=to_features(data,lower,upper)
    f=wcs.wcs_pix2world(ii.T,0)
    f=f.T
    return f


@support_nddata
def integrate(data, wcs=None, mask=None, unit=None, axis=(0)):
    if mask is not None:
        data=fix_mask(data,mask)
    newdata = np.sum(data, axis=axis)
    mask = np.isnan(newdata)
    return NDData(newdata, uncertainty=None, mask=mask, wcs=wcs, meta=None, unit=unit)




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




