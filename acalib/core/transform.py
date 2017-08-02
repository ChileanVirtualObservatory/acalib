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
    inputCont : acalib.Container
        Container with the images to be scaled.

    mayorAxisTemplate : float
        Axis respect the scale will be performed on all the images of the container.

    Returns
    -------
    result : list
        List with the scaled images.
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
    data : list
        List of (M,N,Z) numpy.ndarray images.
    angle : float
        Rotation reference angle that will be applied to all the images.

    Returns
    -------
    result : tuple
        Tuple with the list of rotated images and the list of rotation angles applied to each one.
    """
    rotatedData = []
    angles = []

    for i in np.arange(len(data)):
        prop = fits_props(data[i])
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
    Performs crop and alignment of a list of data cubes.

    Parameters
    ----------
    data : list
        List of astronomical data cubes (numpy.ndarray).

    angles : list
        List of angles (float) to perform alignment.

    Returns
    -------
    result : list
        List of *aligned* astronomical data cubes (numpy.ndarray).
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
        dxl = int((alignedData[i].shape[0] - minShape[0]) / 2)
        dxr = int((alignedData[i].shape[0] + minShape[0]) / 2)
        dyu = int((alignedData[i].shape[1] - minShape[1]) / 2)
        dyd = int((alignedData[i].shape[1] + minShape[1]) / 2)

        alignedData[i] = alignedData[i][dxl:dxr, dyu:dyd]

    return alignedData


def standarize(data):
    """
    Standarize astronomical data cubes in the 0-1 range.

    Parameters
    ----------
    data : numpy.ndarray or astropy.nddata.NDData or or astropy.nddata.NDData
        Astronomical data cube.

    Returns
    -------
    result : tuple
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
    data : numpy.ndarray or astropy.nddata.NDData or or astropy.nddata.NDData
        Astronomical data cube.
    a : float
        Scale value.
    b : float
        Shift value.

    Returns
    -------
    result : numpy.ndarray or astropy.nddata.NDData
        Unstandarized astronomical cube.
    """
    return a*data+b


def add(data, flux, lower, upper):
    """
    Adds flux to a sub-cube of an astronomical data cube.

    Parameters
    ----------
    data : numpy.ndarray or astropy.nddata.NDData or or astropy.nddata.NDData
        Astronomical data cube.
    flux : numpy.ndarray
        Flux added to the cube.
    lower : float
        Lower bound of the sub-cube to which flux will be added.
    upper : float
        Upper bound of the sub-cube to which flux will be added.
    """

    data_slab, flux_slab = matching_slabs(data, flux, lower, upper)
    data[data_slab] += flux[flux_slab]


def denoise(data, threshold):
    """
    Performs denoising of data cube, thresholding over the threshold value.

    Parameters
    ----------
    data : numpy.ndarray or astropy.nddata.NDData or or astropy.nddata.NDData
        Astronomical data cube.
    threshold : float
        Threshold value used for denoising.

    Returns
    -------
    result : numpy.ndarray
        Denoised (thresholded) astronomical data cube.
    """

    elms = data > threshold
    newdata = np.zeros(data.shape)
    newdata[elms] = data[elms]
    return newdata


def fits_props(img):
    """
    Extracts properties information of the astronomical data cube.

    Identifies astronomical objects applying a threshold, detects the larger one that doens't touch the border and returns its properties, approximating it as a (possibily rotated) ellipse.

    Parameters
    ----------
    img : numpy.ndarray
        Astronomical data cube.

    Returns
    -------
    result : dict
        Dictionary with properties of the larger object: *centroid*, *major*, *minor*, *ratio*, *angle*, *area*, *img*, *clr*, *label*, *orig*.
    """
    flt = threshold_otsu(img)
    otsu = img >= flt
    # NOTE: v Not sure if this feature is necessary, if removed, update the documentation of this function.
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

def blit_add(target, source, location=(0,0), operation=lambda t,s: t+s,
        neutral=0):
    """
    Blits an image into another with the given operation (addition by default), allowing source pixels to be left outside the target.

    Parameters
    ----------
    target : numpy.ndarray
        Target image for the blit.
    source : numpy.ndarray
        Image that will be blitted.
    location : tuple
        (x,y) position where the source will be blitted
    operation: function
        function that receives the actual state of the target the source and returns the target modified.
    neutral: int
        The neutral element value for the pixels outside the source.

    Returns
    -------
    target : numpy.ndarray
        Resulting image after blitting.
    """
    # Find the common ranges of both the target and source.
    str_y = max(0,location[0])
    str_x = max(0,location[1])
    end_y = min(target.shape[0],location[0]+source.shape[0])
    end_x = min(target.shape[1],location[1]+source.shape[1])
    str_y_2 = max(0,-location[0])
    str_x_2 = max(0,-location[1])
    end_y_2 = str_y_2+end_y-str_y
    end_x_2 = str_x_2+end_x-str_x
    # Create final source image of the same size than the target and has it on the right position but also the neutral element value outside.
    final_source = np.ones(target.shape)*neutral
    final_source[str_y:end_y,str_x:end_x] = \
        source[str_y_2:end_y_2,str_x_2:end_x_2]
    # Apply the operation
    target = operation(target,final_source)
    return target
