from astropy.nddata import support_nddata, NDData
from acalib import core
import numpy as np

from acalib.upi.axes import opening, features, axes_units

@support_nddata
def noise_level(data,mask=None,unit=None):
    """Compute the RMS of data.
    """
    #TODO: check photutils background estimation for using that if possible
    if unit is None:
        return core.rms(data,mask)
    else:
        return core.rms(data,mask)*unit

@support_nddata
def standarize(data, wcs=None, unit=None, mask=None, meta=None):
    """ Standarize data: data = a * res + b"""
    if mask is not None:
        data = core.fix_mask(data, mask)
    (res, a, b) = core.standarize(data)
    res = NDData(res, uncertainty=None, mask=mask, wcs=wcs, meta=meta, unit=unit)
    return (res, a, b)


@support_nddata
def unstandarize(data, a, b, wcs=None, unit=None, mask=None, meta=None):
    """ unstandarize data: res = a * data + b"""
    if mask is not None:
        data = core.fix_mask(data, mask)
    res = core.unstandarize(data, a, b)
    res = NDData(res, uncertainty=None, mask=mask, wcs=wcs, meta=meta, unit=unit)
    return res


@support_nddata
def add(data, flux, lower=None, upper=None,wcs=None,unit=None,meta=None,mask=None):
    """ Create a new data with the new flux added. 

    Lower and upper are bounds for data. This operation is border-safe and creates a new object at each call.
    Please use the OO version data.add(flux) for modifying the data itself.
    """
    res = data.copy()
    core.add(res, flux, lower, upper)
    return NDData(res, uncertainty=None, mask=mask, wcs=wcs, meta=None, unit=unit)


@support_nddata
def denoise(data, wcs=None, mask=None, unit=None, threshold=0.0):
    """ Simple denoising given a threshold (creates a new object) """
    newdata = core.denoise(data, threshold)
    return NDData(newdata, uncertainty=None, mask=mask, wcs=wcs, meta=None, unit=unit)


@support_nddata
def world_gaussian(data, mu, P, peak, cutoff, wcs=None):
    """ Creates a gaussian flux at mu position (WCS), with P shape, with a maximum value equal to peak,
    and with compact support up to the cutoff contour """
    Sigma = np.linalg.inv(P)
    window = np.sqrt(2 * np.log(peak / cutoff) * np.diag(Sigma)) * axes_units(data, wcs=wcs)
    lower, upper = opening(data, mu, window, wcs=wcs)
    if np.any(np.array(upper - lower) <= 0):
        return None, lower, upper
    feat = features(data, wcs=wcs, lower=lower, upper=upper)
    feat = np.array(feat.columns.values())
    mu = np.array([x.value for x in mu])
    res = core.gaussian_function(mu, P, feat, peak)
    # TODO Not generic
    res = res.reshape(upper[0] - lower[0], upper[1] - lower[1], upper[2] - lower[2])
    return res, lower, upper
