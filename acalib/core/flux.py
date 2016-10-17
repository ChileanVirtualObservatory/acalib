

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
def add_flux(data,flux,lower=None,upper=None):
    """ Adds flux to data. 

    Lower and upper are bounds for data. This operation is border-safe. 
    """
    #if data.ndim!=flux.ndim:
    #    log.error("")

    data_slab,flux_slab=matching_slabs(data,flux,lower,upper)
    data[data_slab]+=flux[flux_slab]

@support_nddata
def denoise(data,wcs=None,mask=None,unit=None,threshold=0.0):
      elms=data>threshold
      newdata=np.zeros(data.shape)
      newdata[elms]=data[elms]
      return NDData(newdata, uncertainty=None, mask=mask,wcs=wcs, meta=None, unit=unit)


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


