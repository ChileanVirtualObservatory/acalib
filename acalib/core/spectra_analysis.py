import numpy as np

# TODO: Not Generic
from astropy import log

from acalib.core import image_analysis
from numpy import NaN, Inf, arange, isscalar, asarray, array

def peakdet(v, delta, x = None):
    """
    Converted from MATLAB script at http://billauer.co.il/peakdet.html
    Returns two arrays
    function [maxtab, mintab]=peakdet(v, delta, x)
    %PEAKDET Detect peaks in a vector
    %        [MAXTAB, MINTAB] = PEAKDET(V, DELTA) finds the local
    %        maxima and minima ("peaks") in the vector V.
    %        MAXTAB and MINTAB consists of two columns. Column 1
    %        contains indices in V, and column 2 the found values.
    %
    %        With [MAXTAB, MINTAB] = PEAKDET(V, DELTA, X) the indices
    %        in MAXTAB and MINTAB are replaced with the corresponding
    %        X-values.
    %
    %        A point is considered a maximum peak if it has the maximal
    %        value, and was preceded (to the left) by a value lower by
    %        DELTA.
    % Eli Billauer, 3.4.05 (Explicitly not copyrighted).
    % This function is released to the public domain; Any use is allowed.
    """
    maxtab = []
    mintab = []

    if x is None:
        x = arange(len(v))

    v = asarray(v)

    if len(v) != len(x):
        log.error('Input vectors v and x must have same length')
        return


    if not isscalar(delta):
        log.error('Input argument delta must be a scalar')
        return

    if delta <= 0:
        log.error('Input argument delta must be positive')
        return

    mn, mx = Inf, -Inf
    mnpos, mxpos = NaN, NaN

    lookformax = True

    for i in arange(len(v)):
        this = v[i]
        if this > mx:
            mx = this
            mxpos = x[i]
        if this < mn:
            mn = this
            mnpos = x[i]

        if lookformax:
            if this < mx-delta:
                maxtab.append((mxpos, mx))
                mn = this
                mnpos = x[i]
                lookformax = False
        else:
            if this > mn+delta:
                mintab.append((mnpos, mn))
                mx = this
                mxpos = x[i]
                lookformax = True

    return array(maxtab), array(mintab)

def spectra_sketch(data, samples, random_state=None):
    """
    Create the sketch spectra using pixel samples.

    Parameters
    ----------
    data : (M,N,Z) numpy.ndarray or astropy.nddata.NDData or astropy.nddata.NDDataRef
        Astronomical data cube.

    samples : Number of pixel samples(int) used for the sketch.

    random_state : (default=None)

    Returns
    -------

     spectra (array) and slices (list)

    """
    # Specific for a FREQ,DEC,RA order
    if random_state is not None:
        random = np.random.RandomState(random_state)
    else:
        random = np.random

    dims = data.shape
    P_x = dims[2]
    P_x_range = range(P_x)
    P_y = dims[1]
    P_y_range = range(P_y)
    frec = dims[0]

    spectra = np.zeros(frec)
    x_ = random.choice(P_x_range, samples, replace=True)
    y_ = random.choice(P_y_range, samples, replace=True)

    pixels = data[:, y_, x_].T
    for pixel in pixels:
        pixel_masked = image_analysis.pixel_processing(pixel)
        spectra += pixel_masked

    spectra = image_analysis.pixel_processing(spectra)

    slices = []
    min_slice = -1
    max_slice = -1
    for i in range(frec - 1):
        if spectra[i] != 0:
            if min_slice == -1:
                min_slice = i
            else:
                if spectra[i + 1] == 0:
                    max_slice = i + 1
                    slices.append(slice(min_slice, max_slice))
                    min_slice = -1
                else:
                    if i == frec - 2:
                        max_slice = i + 1
                        slices.append(slice(min_slice, max_slice))
    return spectra, slices