import numpy as np
from astropy import log
from astropy.nddata import support_nddata, NDData
try:
    from skimage.filters import threshold_adaptive
except:
    from skimage.filter import threshold_adaptive
from skimage.measure import label,regionprops
from skimage.morphology import binary_opening, disk
from skimage.segmentation import clear_border

from ._morph import *

from astropy.table import Table

from .utils import fix_mask, slab
from acalib.core import *

def rms(data, mask=None):
    """Compute the RMS of data. If mask != None, then 
       we use that mask.
    """
    # TODO: check photutils background estimation for using that if possible
    if mask is not None:
        data = fix_mask(data, mask)
    mm = data * data
    rms = np.sqrt(mm.sum() * 1.0 / mm.size)
    return rms

def snr_estimation(data, mask=None, noise=None, points=1000, full_output=False):
    """Heurustic that uses the inflexion point of the thresholded RMS to estimate 
       where signal is dominant w.r.t. noise
    """
    data=fix_mask(data,mask)
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
    """ Returns a numpy array with the integration results. """
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
    
    :param samples: Number of pixel samples used for the sketch.
    :type samples: int
    :returns: ( spectra (array), slices  (list)).
    """
    # Specific for a FREQ,DEC,RA order
    if random_state is not None:
        np.random.seed(random_state)

    dims = data.shape
    P_x = dims[2]
    P_x_range = range(P_x)
    P_y = dims[1]
    P_y_range = range(P_y)
    frec = dims[0]

    spectra = np.zeros(frec)

    for i in range(samples):
        x_ = np.random.choice(P_x_range, 1)
        y_ = np.random.choice(P_y_range, 1)
        pixel = data[:, y_, x_]
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
    acum = _accumulating(pixels)
    diff = _differenting(acum)
    boxing = _segmenting(diff)
    boxing = _erosing(boxing)
    return _masking(boxing, pixels)


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

#Remove NDData support!
@support_nddata
def gaussian_mix(data, prob=0.05, precision=0.02, wcs=None):
    """
    Using a mixture of gaussians make an multiscale segmentation to get the region of interest of a 2D astronomical image.
    
    :param image: Velocity collapsed image
    :returns: list of skimage.measure.regionprops Objects, with detected regions properties
    """
    #TODO: check for wcs != None
    if len(data.shape) > 2:
        log.error("Only 2D images supported")
        raise ValueError("Only 2D images supported")

    image_list = []

    image = data
    image[np.isnan(image)] = 0

    prob = prob
    dims = image.shape
    rows = dims[0]
    cols = dims[1]
    size = np.min([rows, cols])
    precision = size * precision

    image = image.astype('float64')

    w_max = _optimal_w(image, prob)
    diff = (image - np.min(image)) / (np.max(image) - np.min(image))

    tt = w_max * w_max
    if tt % 2 == 0:
        tt += 1
    g = threshold_adaptive(diff, tt, method='mean', offset=0)

    r = w_max / 2
    rMin = 2 * np.round(precision)

    while (r > rMin):
        background = np.zeros((rows, cols))
        selem = disk(r)
        sub = binary_opening(g, selem)
        sub = clear_border(sub)
        sub = label(sub)
        fts = regionprops(sub)

        #image_list.append(NDData(sub, wcs=wcs))
        # Non NNData version (without wcs... lets check if it pass)
        image_list.append(sub)

        if len(fts) > 0:
            for props in fts:
                C_x, C_y = props.centroid

                radius = props.equivalent_diameter / 2.
                kern = 0.01 * np.ones((2 * radius, 2 * radius))
                krn = _kernelsmooth(x=np.ones((2 * radius, 2 * radius)), kern=kern)
                krn = np.exp(np.exp(krn))
                if np.max(krn) > 0:
                    krn = (krn - np.min(krn)) / (np.max(krn) - np.min(krn))
                    background = _kernel_shift(background, krn, C_x, C_y)
        if np.max(background) > 0:
            background = (background - np.min(background)) / (np.max(background) - np.min(background))
            diff = diff - background
        diff = (diff - np.min(diff)) / (np.max(diff) - np.min(diff))
        tt = r * r
        if tt % 2 == 0:
            tt += 1
        g = threshold_adaptive(diff, tt, method='mean', offset=0)
        r = np.round(r / 2.)

    return image_list


def _optimal_w(image, p=0.05):
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
        tt = radius * radius
        if tt % 2 == 0:
            tt += 1

        g = threshold_adaptive(f, tt, method='mean', offset=0)
        ov = _bg_fg(f, g, bg, fg)
        if (ov < min_ov):
            w = radius
            min_ov = ov

        radius += inc
    return w


def _bg_fg(f, g, bg, fg):
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
def vel_stacking(data, data_slice):
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
    # wcs = wcs.dropaxis(2)

    return stacked  #NDData(stacked, uncertainty=uncertainty, mask=mask,wcs=wcs, meta=meta, unit=unit)


### DEPRECATED ####


# def ndslice(ndd, lower, upper):
#    """
#    N-Dimensional slicing.
#
#    Arguments:
#        ndd   -- an astropy.nddata.NDDataArray object.
#        lower -- n-dimensional point as an n-tuple.
#        upper -- n-dimensional point as an n-tuple.
#
#    Returns:
#        A sliced astropy.nddata.NDDataArray object.
#
#    """
#    lower = lower if lower is not None else np.zeros(ndd.ndim)
#    upper = upper if upper is not None else ndd.shape
#    return ndd[[slice(min(a,b), max(a,b)+1) for a,b in zip(lower, upper)]]
#
# def adjust_index(relative, origin):
#    """
#    Adjusts an index relative to a subarray to an absolute
#    index in the superarray.
#
#    Arguments:
#        origin   -- an n-dimensional index of a point as an n-tuple.
#                    It should be the origin from which the relative
#                    index was computed.
#        relative -- an n-dimensional index of a point as an n-tuple.
#                    The index to be adjusted.
#
#    Returns:
#        The relative index adjusted to the superarray as an n-tuple.
#    """
#    return tuple(np.array(origin) + np.array(relative))

# def index_of_max(ndd, lower=None, upper=None):
#    """
#    Index of maximum value in an m-dimensional subarray from
#    an n-dimensional array, specified by lower and upper.
#
#    Arguments:
#        ndd   -- an astropy.nddata.NDDataArray object.
#        lower -- n-dimensional point as an n-tuple.
#        upper -- n-dimensional point as an n-tuple.
#
#    Returns:
#        A tuple with the maximum value found in the m-dimensional
#        subarray and its index in the n-dimensional superarray.
#
#    """
#    ndd = ndslice(ndd, lower, upper)
#    index = np.unravel_index(ndd.data.argmax(), ndd.data.shape)
#    value = ndd.data[index]
#    return (value, adjust_index(index, lower))

# def index_of_min(ndd, lower=None, upper=None):
#    """
#    Index of minimum value in an m-dimensional subarray from
#    an n-dimensional array, specified by lower and upper.
#
#    Arguments:
#        ndd   -- an astropy.nddata.NDDataArray object.
#        lower -- n-dimensional point as an n-tuple.
#        upper -- n-dimensional point as an n-tuple.
#
#    Returns:
#        A tuple with the minimum value found in the m-dimensional
#        subarray and its index in the n-dimensional superarray.
#
#    """
#    ndd = ndslice(ndd, lower, upper)
#    index = np.unravel_index(ndd.data.argmin(), ndd.data.shape)
#    value = ndd.data[index]
#    return (value, adjust_index(index, lower))
