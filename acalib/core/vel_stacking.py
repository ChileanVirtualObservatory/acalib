from astropy.nddata import support_nddata, NDData
import numpy as np
from astropy import log

@support_nddata
def vel_stacking(data,data_slice, wcs=None, mask=None,uncertainty=None, meta=None, unit=None):
    """
        Create an image collapsing the frecuency axis
        
        :param data_slice: Sector to be collapsed
        :type data_slice: slice 
        :returns: image (NDData): 2D-Array with the stacked cube.

    """
    if len(data.shape) != 3:
        log.error("Cube needs to be a 3D array")        
        raise ValueError("Cube needs to be a 3D array")

    dims = data.shape
    subcube = data[data_slice, :,:]
    stacked = np.sum(subcube,axis=0)
    wcs = wcs.dropaxis(3)

    return NDData(stacked, uncertainty=uncertainty, mask=mask,wcs=wcs, meta=meta, unit=unit)