
def gaussian_function(mu,P,feat,peak):
    """ Generates an n-dimensional Gaussian using the feature matrix feat,
    centered at mu, with precision matrix P and with intensity peak.
    """
    cent_feat=np.empty_like(feat)
    for i in range(len(mu)):
       cent_feat[i]=feat[i] - mu[i]
    qform=(P.dot(cent_feat))*cent_feat
    quad=qform.sum(axis=0)
    res=np.exp(-quad/2.0)
    res=peak*(res/res.max())
    return res

@support_nddata
def world_gaussian(data,wcs,mu,P,peak,cutoff):
   """ Creates a gaussian flux at mu position (WCS), with P shape, with a maximum value equal to peak, 
   and with compact support up to the cutoff contour """
   Sigma=np.linalg.inv(P)
   window=np.sqrt(2*np.log(peak/cutoff)*np.diag(Sigma))
   lower,upper=fov_to_index(data,wcs,mu,window)
   if np.any(np.array(upper-lower)<=0):
       return None,lower,upper
   feat=world_features(data,wcs,lower,upper)
   res=gaussian_function(mu,P,feat,peak)
   # TODO Not generic
   res=res.reshape(upper[0]-lower[0],upper[1]-lower[1],upper[2]-lower[2])
   return res,lower,upper


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

