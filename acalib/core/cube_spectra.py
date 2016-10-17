from astropy.nddata import support_nddata, NDData
from astropy import log 
import numpy as np

from acalib.core._morph import *

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
