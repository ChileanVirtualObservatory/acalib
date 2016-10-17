from astropy.nddata import support_nddata, NDData
from astropy import log 
from acalib.core._morph import *
import numpy as np


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
    wcs = wcs.dropaxis(2)

    return NDData(stacked, uncertainty=uncertainty, mask=mask,wcs=wcs, meta=meta, unit=unit)


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


@support_nddata
def integrate(data, wcs=None, mask=None, unit=None, axis=(0)):
    if mask is not None:
        data=fix_mask(data,mask)
    newdata = np.sum(data, axis=axis)
    mask = np.isnan(newdata)
    return NDData(newdata, uncertainty=None, mask=mask, wcs=wcs, meta=None, unit=unit)




@support_nddata
def cube_spectra(data,samples, random_state = None):
    """
    Create the spectra using pixel samples.
    
    :param samples: Number of pixel samples used for the sketch.
    :type samples: int
    :returns: ( spectra (array), slices  (list)).
    """

    if random_state is not None:
        np.random.seed(random_state)

    dims = data.shape
    P_x = dims[2]
    P_x_range = range(P_x)
    P_y = dims[1]
    P_y_range = range(P_y)
    frec = dims[0]

    spectra = np.zeros(frec)

    for i in xrange(samples):
        x_ = np.random.choice(P_x_range,1)
        y_ = np.random.choice(P_y_range,1)
        pixel = data[:, y_, x_] 
        pixel_masked = _pixel_processing(pixel)
        spectra += pixel_masked
    spectra = _pixel_processing(spectra)

    slices = []
    min_slice = -1
    max_slice = -1
    for i in range(frec-1):
        if spectra[i] != 0:
            if min_slice == -1:
                min_slice = i
            else:
                if spectra[i+1] == 0:
                    max_slice = i+1
                    slices.append(slice(min_slice,max_slice))
                    min_slice = -1
                else:
                    if i == frec-2:
                        max_slice = i+1
                        slices.append(slice(min_slice,max_slice))

    return spectra,slices

def _pixel_processing(pixels):
    pixels = pixels.astype(np.float64)
    acum = _accumulating(pixels)
    diff = _differenting(acum)
    boxing = _segmenting(diff)
    boxing = _erosing(boxing)
    return _masking(boxing,pixels)


def _accumulating(pixels):
    return np.cumsum(pixels)

def _differenting(cumPixels):
    d = diff(cumPixels)
    return d

def _segmenting(diff):
    
    boxing = seg(diff)
    
    return boxing

def _erosing(boxing):
    
    boxing = eros(boxing)

    return boxing

def _masking(boxing, pixels):
    return boxing.reshape(-1)*pixels.reshape(-1)
