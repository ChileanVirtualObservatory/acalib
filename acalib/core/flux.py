from astropy.nddata import support_nddata, NDData
import numpy as np

def fix_mask(data,mask):
    ismasked=isinstance(data,np.ma.MaskedArray)
    if ismasked and mask is None:
        return data
    else:
       return np.ma.MaskedArray(data,mask)

def _standarize(data):
    y_min=data.min()
    res=data - y_min
    y_fact=res.sum()
    res=res/y_fact
    return (res,y_fact,y_min)

def _unstandarize(data,a,b):
    return data*a + b

def _add(data,flux,lower,upper):
    data_slab,flux_slab=matching_slabs(data,flux,lower,upper)
    data[data_slab]+=flux[flux_slab]

def _denoise(data,threshold):
    elms=data>threshold
    data[elms]=data[elms]

@support_nddata
def standarize(data,wcs=None,unit=None,mask=None,meta=None):
    """ Standarize data: data = a * res + b"""
    if mask is not None:
        data=fix_mask(data,mask)
    (res,a,b)=_standarize(data)
    res=NDData(res, uncertainty=None, mask=mask,wcs=wcs, meta=meta, unit=unit)
    return (res,a,b)

@support_nddata
def unstandarize(data,wcs=None,unit=None,mask=None,meta=None):
    """ unstandarize data: res = a * data + b"""
    if mask is not None:
        data=fix_mask(data,mask)
    res=_unstandarize(data,a,b)
    res=NDData(res, uncertainty=None, mask=mask,wcs=wcs, meta=meta, unit=unit)
    return res

@support_nddata
def add(data,flux,lower=None,upper=None):
    """ Create a new data with the new flux added. 

    Lower and upper are bounds for data. This operation is border-safe and creates a new object at each call.
    Please use the OO version data.add(flux) for modifying the data itself.
    """
    res=data.copy()
    _add(res,flux,lower,upper)
    return res

@support_nddata
def denoise(data,wcs=None,mask=None,unit=None,threshold=0.0):
      """ Simple denoising given a threshold (creates a new object) """
      newdata=data.copy()
      _denoise(newdata,threshold)
      return NDData(newdata, uncertainty=None, mask=mask,wcs=wcs, meta=None, unit=unit)


