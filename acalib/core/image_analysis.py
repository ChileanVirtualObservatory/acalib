import numpy as np
from astropy import log
from skimage.filters import threshold_local
from skimage.measure import label,regionprops
from .utils import fix_mask, slab

from ._morph import differenceImpl, segmentationImpl, erosionImpl

from acalib.core import *

def fits_props(img):
    """
    Extracts properties information of the astronomical data cube.

    Parameters
    ----------
    img : numpy.ndarray
        Astronomical data cube.

    Returns
    -------
    result : dict
        Dictionary with properties of the image: *centroid*, *major*, *minor*, *ratio*, *angle*, *area*, *img*, *clr*, *label*, *orig*.
    """
    flt = threshold_otsu(img)
    otsu = img >= flt
    clr = clear_border(otsu)

    # label image regions
    label_image, nlabel = label(clr, return_num=True)
    borders = np.logical_xor(otsu, clr)
    label_image[borders] = -1

    props = regionprops(label_image)

    ratios = []
    areas = []

    for i in props:
        if i.major_axis_length > 0:
            ratios.append(i.minor_axis_length / i.major_axis_length)
            areas.append(i.area)
        else:
            ratios.append(0)
            areas.append(0)

    if len(props) > 1:
        pos = np.argmax(areas)
    else:
        pos = 0

    properties = {'centroid': props[pos].centroid, 'major': props[pos].major_axis_length,
                  'minor': props[pos].minor_axis_length, 'ratio': ratios[pos],
                  'angle': props[pos].orientation, 'area': props[pos].area, 'img': props[pos].image,
                  'clr': clr, 'label': label_image, 'orig': img}
    return properties


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

def pixel_processing(pixels):
    pixels = pixels.astype(np.float64)
    acum = np.cumsum(pixels)
    diff = differenceImpl(acum)
    boxing = segmentationImpl(diff)
    boxing = erosionImpl(boxing)
    return _masking(boxing,pixels)

def _masking(boxing, pixels):
    return boxing.reshape(-1) * pixels.reshape(-1)

def optimal_w(image, p=0.05):
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


def kernelsmooth(x, kern, norm=True):
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

def kernel_shift(back, kernel, x, y):
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
