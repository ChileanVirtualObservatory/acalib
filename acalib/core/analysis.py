import numpy as np
from astropy import log
from astropy.nddata import support_nddata, NDDataRef
from skimage.filters import threshold_local
from skimage.measure import label,regionprops

from ._morph import differenceImpl, segmentationImpl, erosionImpl

from astropy.table import Table

from .utils import fix_mask, slab
from acalib.core import *

def rms(data, mask=None):
    """
    Compute the RMS of data. If mask != None, then we use that mask.

    Parameters
    ----------
    data : (M,N,Z) numpy.ndarray or astropy.nddata.NDData or or astropy.nddata.NDDataRef
        Astronomical data cube.

    mask : numpy.ndarray (default = None)

    Returns
    -------
    RMS of the data (float)
    """
    # TODO: check photutils background estimation for using that if possible
    if mask is not None:
        data = fix_mask(data, mask)
    mm = data * data
    rms = np.sqrt(mm.sum() * 1.0 / mm.size)
    return rms

def snr_estimation(data, mask=None, noise=None, points=1000, full_output=False):
    """
    Heurustic that uses the inflexion point of the thresholded RMS to estimate where signal is dominant w.r.t. noise

    Parameters
    ----------
    data : (M,N,Z) numpy.ndarray or astropy.nddata.NDData or astropy.nddata.NDDataRef
        Astronomical data cube.

    mask : numpy.ndarray (default = None)

    noise : float (default=None)
        Noise level, if not given will use rms of the data.

    points : (default=1000)

    full_output : boolean (default=False)
        Gives verbose results if True

    Returns
    --------

    "Signal to Noise Radio" value

    """
    if noise is None:
        noise = rms(data, mask)
    x = []
    y = []
    n = []
    sdata = data[data > noise]
    for i in range(1, int(points)):
        val = 1.0 + 2.0 * i / points
        sdata = sdata[sdata > val * noise]
        if sdata.size < 2:
            break
        n.append(sdata.size)
        yval = sdata.mean() / noise
        x.append(val)
        y.append(yval)
    y = np.array(y)
    v = y[1:] - y[0:-1]
    p = v.argmax() + 1
    snrlimit = x[p]
    if full_output == True:
        return snrlimit, noise, x, y, v, n, p
    return snrlimit


def integrate(data, mask=None, axis=(0)):
    """
    Sums the slices of a cube of data given an axis.

    Parameters
    ----------
    data : (M,N,Z) numpy.ndarray or astropy.nddata.NDData or astropy.nddata.NDDataRef
        Astronomical data cube.

    mask : numpy.ndarray (default = None)

    axis : int (default=(0))

    Returns
    -------
     A numpy array with the integration results.

    """
    if mask is not None:
        data = fix_mask(data, mask)
    newdata = np.sum(data, axis=axis)
    mask = np.isnan(newdata)
    return newdata


def get_shape(data, intensity_image, wcs=None):
    objs_properties = []
    fts = regionprops(data, intensity_image=intensity_image)

    for obj in fts:
        if wcs:
            matrix = wcs.pixel_scale_matrix
            deg_per_pix_x = matrix[0, 0]
            deg_per_pix_y = matrix[1, 1]

        centroid = wcs.celestial.all_pix2world(obj.centroid[0], obj.centroid[1], 1) if wcs else obj.centroid
        centroid_ra = centroid[0] if wcs else centroid[0]
        centroid_dec = centroid[1] if wcs else centroid[1]
        major_axis = abs(deg_per_pix_x) * obj.major_axis_length if wcs else obj.major_axis_length
        minor_axis = abs(deg_per_pix_x) * obj.minor_axis_length if wcs else obj.minor_axis_length
        area = obj.area * abs(deg_per_pix_x) if wcs else obj.area
        eccentricity = obj.eccentricity
        solidity = obj.solidity
        filled = obj.area / obj.filled_area

        objs_properties.append((centroid_ra, centroid_dec, major_axis, minor_axis, area,
                                eccentricity, solidity, filled, obj.max_intensity, obj.min_intensity,
                                obj.mean_intensity))

    return objs_properties


#TODO: try not to use nddata and Table
@support_nddata
def measure_shape(data, labeled_images, min_freq=None, max_freq=None, wcs=None):
    """ Measure a few statistics from labeled images """
    # TODO: Document this function
    objects = list()
    intensity_image = data
    for image in labeled_images:
        objs_properties = get_shape(image, intensity_image)
        objects.extend(objs_properties)

    if len(objects) == 0:
        return Table()

    names = ["CentroidRa", "CentroidDec", "MajorAxisLength", "MinorAxisLength",
             "Area", "Eccentricity", "Solidity", "FilledPercentaje", "MaxIntensity", "MinIntensity", "AverageIntensity"]

    meta = {"name": "Object Shapes"}

    if min_freq is not None:
        meta["min_freq_hz"] = min_freq

    if max_freq is not None:
        meta["max_freq_hz"] = max_freq

    t = Table(rows=objects, names=names, meta=meta)
    return t


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
        pixel_masked = _pixel_processing(pixel)
        spectra += pixel_masked
    
    spectra = _pixel_processing(spectra)

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

def _pixel_processing(pixels):
    pixels = pixels.astype(np.float64)
    acum = np.cumsum(pixels)
    diff = differenceImpl(acum)
    boxing = segmentationImpl(diff)
    boxing = erosionImpl(boxing)
    return _masking(boxing,pixels)

def _masking(boxing, pixels):
    return boxing.reshape(-1) * pixels.reshape(-1)

def index_mesh(data, lower=None, upper=None):
    """ Create an meshgrid from indices """
    sl = slab(data, lower, upper)
    dim = data.ndim
    slices = []
    for i in range(dim):
        slices.append(slice(sl[i].start, sl[i].stop))
    retval = np.mgrid[slices]
    return retval


def index_features(data, lower=None, upper=None):
    """ Creates an array with indices in features format """
    msh = index_mesh(data, lower, upper)
    dim = data.ndim
    ii = np.empty((dim, int(msh.size / dim)))
    for i in range(dim):
        ii[dim - i - 1] = msh[i].ravel()
    return ii

def _optimal_w(image, p=0.05):
    # Calculate the optimal window size for the image segmentation given a quantile.
    # It expand the radious until it reaches the best segmentation.


    # radiusMin, radius Max and inc in percentages of the image size, p as [0,1] value, image is the original version
    radiusMin = 5
    radiusMax = 40
    inc = 1

    f = (image - np.min(image)) / (np.max(image) - np.min(image))
    dims = f.shape
    rows = dims[0]
    cols = dims[1]

    maxsize = np.max([rows, cols])
    imagesize = cols * rows
    radius_thresh = np.round(np.min([rows, cols]) / 4.)
    unit = np.round(maxsize / 100.)

    radiusMin = radiusMin * unit
    radiusMax = radiusMax * unit
    radiusMax = int(np.min([radiusMax, radius_thresh]))
    radius = radiusMin
    inc = inc * unit

    bg = np.percentile(f, p * 100)
    fg = np.percentile(f, (1 - p) * 100)
    min_ov = imagesize

    while (radius <= radiusMax):
        tt = int(radius * radius)
        if tt % 2 == 0:
            tt += 1

        adaptive_threshold = threshold_local(f, tt, method='mean', offset=0)#(f, tt, offset=0)
        g = f > adaptive_threshold

        ov = _bg_fg(f, g, bg, fg)
        if (ov < min_ov):
            w = radius
            min_ov = ov

        radius += inc
    return w


def _bg_fg(f, g, bg, fg):
    # Calculate the backgorund and foreground distribution
    #
    dims = f.shape
    rows = dims[0]
    cols = dims[1]
    fp = 0
    fn = 0
    for rowID in range(rows):
        for colId in range(cols):
            if g[rowID][colId] == True:
                if (np.abs(f[rowID][colId] - bg) < np.abs(f[rowID][colId] - fg)):
                    fp += 1
            elif g[rowID][colId] == False:
                if (np.abs(f[rowID][colId] - bg) > np.abs(f[rowID][colId] - fg)):
                    fn += 1
    overall = fp + fn
    return overall


def _kernelsmooth(x, kern, norm=True):
    # how many rows/cols of zeroes are used to pad.
    width = kern.shape[0]
    pad = int(width / 2.)

    # record the width and height the input data matrix
    x_w = x.shape[0]
    x_h = x.shape[1]

    if norm:
        k = kern / np.sum(abs(kern))
    else:
        k = kern

    # Padding with zeros
    x_pad = np.lib.pad(x, ((pad, pad), (pad, pad)), 'constant')

    # Pre-allocate the final (smoothed) data matrix
    s = np.zeros((x_h, x_w))

    # Pre-allocate a temporary matrix for the iterative calculations
    temp = np.zeros((width, width))

    # Loop through the data to apply the kernel.
    for col in range(x_w):
        for row in range(x_h):
            temp = x_pad[row:(row + width), col:(col + width)]
            s[row][col] = np.sum(k * temp)

    return s


def _kernel_shift(back, kernel, x, y):
    rows_back = back.shape[0]
    cols_back = back.shape[1]
    rowsKernel = kernel.shape[0]
    colsKernel = kernel.shape[1]
    rowInit = int(x - rowsKernel / 2)
    colInit = int(y - colsKernel / 2)

    for row in range(rowsKernel):
        for col in range(colsKernel):
            if (rowInit + row < rows_back - 1) and (colInit + col < cols_back - 1):
                back[rowInit + row][colInit + col] = kernel[row][col]

    return back


# TODO: This is non-generic, uses the axis=0!
@support_nddata
def vel_stacking(data,data_slice,wcs=None,uncertainty=None, mask=None, meta=None, unit=None):
    """
    Create an image collapsing the frecuency axis

    Parameters
    ----------
    data : numpy.ndarray or astropy.nddata.NDData or astropy.nddata.NDDataRef
        Astronomical 2D image

    slice : slice object
        Sector to be collapsed

    Returns
    -------
    image (NDDataRef): 2D-Array with the stacked cube.

    """
    if len(data.shape) != 3:
        log.error("Cube needs to be a 3D array")
        raise ValueError("Cube needs to be a 3D array")
    dims = data.shape
    subcube = data[data_slice, :,:]
    stacked = np.sum(subcube,axis=0)
    if wcs:
        wcs = wcs.dropaxis(2)

        return NDDataRef(stacked, uncertainty=uncertainty, mask=mask,wcs=wcs, meta=meta, unit=unit)
    else:
        return stacked