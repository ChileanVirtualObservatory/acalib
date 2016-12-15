import scipy.ndimage as scnd
import numpy as np


from skimage.filters import threshold_otsu
from skimage.measure import label
from skimage.segmentation import clear_border
from skimage.measure import regionprops


from . import utils


def scale(inputCont, majorAxisTemplate):
    """
    Performs an scale of the images in the container acording to the indicated mayor axis.

    Parameters
    ----------
    inputCont: acalib.Container
        Container with the images to be scaled.

    mayorAxisTemplate: float
        Axis respect the scale will be performed on all the images of the container.

    Returns
    -------
    result: Python list with the scaled images. 
    """
    scaledData = []

    for i in np.arange(len(inputCont.images)):
        prop = fits_props(inputCont.images[i].data)
        sc = majorAxisTemplate / prop['major']
        scaledData.append(scnd.zoom(prop['orig'], sc))
    return scaledData


def rotate(data, angle):
    """
    Performs an angle rotation over a list of images
    

    Parameters
    ----------
    data: list
        List of (M,N,Z) numpy.ndarray images.
    angle: float
        Rotation reference angle that will be applied to all the images.

    Returns
    -------
    result: tuple
        Tuple with the list of rotated images and the list of rotation angles applied to each one.  
    """
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
    """
    Parameters
    ----------
    data: list
        

    angles: list


    Returns
    -------
    result: 
    """
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
    """
    Standarize astronomical data cubes in the 0-1 range.

    Parameters
    ----------
    data: numpy.ndarray or astropy.nddata.NDData
        Astronomical data cube.

    Returns
    -------
    result: tuple
        Tuple containing the standarized numpy.ndarray or astropy.nddata.NDData cube, the factor scale y_fact and the shift y_min.
    """
    y_min = data.min()
    res = data - y_min
    y_fact = res.sum()
    res = res / y_fact
    return (res, y_fact, y_min)


def unstandarize(data, a, b):
    """
    Unstandarize the astronomical data cube: :math:`a \cdot data + b`.

    Parameters
    ----------
    data: numpy.ndarray or astropy.nddata.NDData
        Astronomical data cube.
    a: float
        Scale value.
    b: float
        Shift value.

    Returns
    ------- 
    result: numpy.ndarray or astropy.nddata.NDData
        Unstandarized astronomical cube.
    """
    return a*data+b


def add(data, flux, lower, upper):
    """

    """
    data_slab, flux_slab = matching_slabs(data, flux, lower, upper)
    data[data_slab] += flux[flux_slab]

def denoise(data, threshold):
    """
    Parameters
    ----------

    Returns
    -------
    """
    elms = data > threshold
    newdata = np.zeros(data.shape)
    newdata[elms] = data[elms]
    return newdata


def fits_props(img):
    """

    Parameters
    ----------

    Returns
    -------
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
        ratios.append(i.minor_axis_length / i.major_axis_length)
        areas.append(i.area)

    if len(props) > 1:
        pos = areas.index(max(areas))
    else:
        pos = 0

    properties = {'centroid': props[pos].centroid, 'major': props[pos].major_axis_length,
                  'minor': props[pos].minor_axis_length, 'ratio': ratios[pos],
                  'angle': props[pos].orientation, 'area': props[pos].area, 'img': props[pos].image,
                  'clr': clr, 'label': label_image, 'orig': img}

    return properties


def img_props(img):
    """

    Parameters
    ----------

    Returns
    -------
    result:
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
        ratios.append(i.minor_axis_length / i.major_axis_length)
        areas.append(i.area)

    if len(props) > 1:
        pos = areas.index(max(areas))
    else:
        pos = 0

    properties = {'centroid': props[pos].centroid, 'major': props[pos].major_axis_length,
                  'minor': props[pos].minor_axis_length, 'ratio': ratios[pos],
                  'angle': props[pos].orientation, 'area': props[pos].area, 'img': props[pos].image,
                  'clr': clr, 'label': label_image, 'orig': img}

    return properties
