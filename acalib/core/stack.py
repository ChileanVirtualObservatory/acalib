import numpy as np
from scipy.stats import signaltonoise
from skimage.filters import threshold_otsu
from skimage.measure import label
from skimage.segmentation import clear_border
from skimage.measure import regionprops

from acalib.core import *


def stacking(template_data, data_cont):
    tprops = fits_props(template_data)
    scaled = scale(data_cont, tprops['major'])
    rotated, angles = rotate(scaled, tprops['angle'])
    aligned = crop_and_align(rotated, angles)
    result = np.mean(aligned, axis=0)
    return result, signaltonoise(result)

def img_props(img):
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
