from astropy.nddata import support_nddata, NDData
from astropy import log 
import numpy as np

@support_nddata
def cube_spectra(data,samples, random_state=None):
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
    acum = _accumulating(pixels)
    diff = _differenting(acum)
    boxing = _segmenting(diff)
    boxing = _erosing(boxing)
    return _masking(boxing,pixels)


def _accumulating(pixels):
    return np.cumsum(pixels)

def _differenting(cumPixels):
    n = len(cumPixels)
    diff = np.zeros(n)
    diff[0] = cumPixels[0]
    for i in xrange(1,n):
        diff[i] = cumPixels[i] - diff[i-1]
    return diff 

def _segmenting(diff):
    n = len(diff)
    boxing = np.zeros(n)
    
    for i in xrange(1,n-1):
        boxing[i] = 1
        if (
            (diff[i] < diff[i-1]) and (diff[i] < diff[i+1])

            or

            (diff[i] > diff[i-1]) and (diff[i] > diff[i+1])
            ):
            boxing[i] = 0

    return boxing

def _erosing(boxing):
    n = len(boxing)
    blocking = np.zeros(n)

    for i in xrange(1,n-1):
        blocking[i] = boxing[i]

        if ( boxing[i-1] == 0 and boxing[i] == 1 and boxing[i+1] == 0 ):
            blocking[i] = 0

    boxing = np.copy(blocking)
    for i in xrange(1, n-1):
        if (blocking[i-1] == 0 and blocking[i]== 1):
            boxing[i-1] = 1
        if (blocking[i]==1 and blocking[i+1]==0):
            boxing[i+1] = 1

    return boxing

def _masking(boxing, pixels):
    n1 = len(boxing)
    n2 = len(pixels)

    if n1 == n2:
        n = n1
        output = np.zeros(n)
        for i in xrange(n):
            output[i] = boxing[i]* pixels[i]
        return output
    else:
        log.error("boxing and pixels has different length")        
        raise ValueError("boxing and pixels has different length")
