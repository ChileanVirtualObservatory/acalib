import scipy.ndimage as scnd
import numpy as np

from .stack import img_props, fits_props
from .utils import matching_slabs


def scale(inputCont, majorAxisTemplate):
    scaledData = []

    for i in np.arange(len(inputCont.images)):
        prop = fits_props(inputCont.images[i].data)
        scale = majorAxisTemplate / prop['major']
        scaledData.append(scnd.zoom(prop['orig'], scale))
    return scaledData


def rotate(data, angle):
    rotatedData = []
    angles = []

    for i in np.arange(len(data)):
        prop = img_props(data[i])
        angles.append(angle - prop['angle'])
        rotatedData.append(scnd.rotate(data[i], angles[-1], reshape=True))
    return rotatedData, angles


def _rotation_limits(img, angle):
    if angle > 0:
        cx, cy = np.nonzero(np.array(img.T))
    else:
        cx, cy = np.nonzero(np.array(img))

    upper = (cx[0], cy[0])
    lower = (cx[-1], cy[-1])

    return upper, lower


def crop_and_align(data, angles):
    alignedData = []
    shapes = []

    for i in np.arange(len(data)):
        upper, lower = _rotation_limits(data[i], angles[i])
        crop = data[i][upper[1]:lower[1], upper[1]:lower[1]]
        shapes.append(list(crop.shape))

        alignedData.append(crop)

    minShape = tuple(np.amin(shapes, axis=0))

    for i in np.arange(len(alignedData)):
        dxl = (alignedData[i].shape[0] - minShape[0]) / 2
        dxr = (alignedData[i].shape[0] + minShape[0]) / 2
        dyu = (alignedData[i].shape[1] - minShape[1]) / 2
        dyd = (alignedData[i].shape[1] + minShape[1]) / 2

        alignedData[i] = alignedData[i][dxl:dxr, dyu:dyd]

    return alignedData


def standarize(data):
    y_min = data.min()
    res = data - y_min
    y_fact = res.sum()
    res = res / y_fact
    return (res, y_fact, y_min)


def unstandarize(data, a, b):
    return data * a + b


def add(data, flux, lower, upper):
    data_slab, flux_slab = matching_slabs(data, flux, lower, upper)
    data[data_slab] += flux[flux_slab]

def denoise(data, threshold):
    elms = data > threshold
    newdata = np.empty_like(data)
    newdata[elms] = data[elms]
    return newdata
