
import numpy as np
import matplotlib.pyplot as plt


def fix_limits(data,vect):
    if isinstance(vect,tuple):
        vect=np.array(vect)
    vect=vect.astype(int)
    low=vect < 0
    up=vect > data.shape
    if vect.any():
        vect[low]=0
    if vect.any():
        vect[up]=np.array(data.shape)[up]
    return vect


def slab(data,lower=None,upper=None):
    if lower is None:
        lower=np.zeros(data.ndim)
    if upper is None:
        upper=data.shape
    lower=fix_limits(data,lower)
    upper=fix_limits(data,upper)
    m_slab=[]
    for i in range(data.ndim):
       m_slab.append(slice(lower[i],upper[i]))
    return m_slab

def matching_slabs(data,flux,lower,upper):
    data_slab=slab(data,lower,upper)
    flow=np.zeros(flux.ndim)
    fup=np.array(flux.shape)
    for i in range(data.ndim):
       if data_slab[i].start == 0:
          flow[i] = flux.shape[i] - data_slab[i].stop
       if data_slab[i].stop == data.shape[i]:
          fup[i] = data_slab[i].stop - data_slab[i].start
    flux_slab=slab(flux,flow,fup)
    return data_slab,flux_slab

def add_flux(data,flux,lower=None,upper=None):
    #if data.ndim!=flux.ndim:
    #    log.error("")
    data_slab,flux_slab=matching_slabs(data,flux,lower,upper)
    data[data_slab]+=flux[flux_slab]


def create_gauss(mu,P,feat,peak):
    cent_feat=np.empty_like(feat)
    for i in range(len(mu)):
       cent_feat[i]=feat[i] - mu[i]
    qform=(P.dot(cent_feat))*cent_feat
    quad=qform.sum(axis=0)
    res=np.exp(-quad/2.0)
    res=peak*(res/res.max())
    return res


def create_mould(P,delta):
    """This function creates a gaussian mould, using the already computed values of 
    """
    # TODO Can we use index_features
    n=len(delta)
    ax=[]
    elms=[]
    for i in range(n):
        lin=np.arange(-delta[i]-0.5,delta[i]+0.5)
        elms.append(len(lin))
        ax.append(lin)
    grid=np.meshgrid(*ax,indexing='ij')
    feat=np.empty((n,np.product(elms)))
    for i in range(n):
        feat[i]=grid[i].ravel()
    mould=create_gauss(np.zeros(n),P,feat,1)
    mould=mould.reshape(*elms)
    return(mould)

def estimate_rms(data):
    mm=data * data
    if isinstance(mm,np.ma.MaskedArray):
        rms=np.sqrt(mm.sum()*1.0/mm.count())
    else:
        rms=np.sqrt(mm.sum()*1.0/mm.size)
    return rms


if __name__ == '__main__':
    # Slab and AddFlux test
    a=np.random.random((20,20,20))
    sl=slab(a,(-5,4,5),(15,25,10))
    print(sl)
    b=100*np.random.random((10,10,10))
    add_flux(a,b,(15,-5,7),(25,5,17))
    c=np.where(a>1.0)
    print(str(c[0].size)+" should be near 250")
    # Mould test
    P=np.array([[0.05,0.01,0],[0.01,0.07,0.03],[0,0.03,0.09]])
    delta=[10,15,20]
    mould=create_mould(P,delta)
    plt.imshow(mould.sum(axis=(0)))
    plt.show()




