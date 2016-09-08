from astropy.nddata import support_nddata, NDData
import numpy as np

@support_nddata
def vel_stacking(data,data_slice, wcs=None, mask=None,uncertainty=None, meta=None, unit=None):
    """
        Create an image collapsing the frecuency axis
        
        :param data_slice: Sector to be collapsed
        :type data_slice: slice 
        :returns: image (NDData): 2D-Array with the stacked cube.

    """
    dims = data.shape
    subcube = data[data_slice, :,:]
    stacked = np.sum(subcube,axis=0)

    return NDData(stacked, uncertainty=uncertainty, mask=mask,wcs=wcs, meta=meta, unit=unit)