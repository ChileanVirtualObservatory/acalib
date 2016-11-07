import numpy as np
from scipy.stats import signaltonoise

from .transform import fits_props, scale, rotate, crop_and_align


def stacking(template_data, data_cont):
    tprops = fits_props(template_data)
    scaled = scale(data_cont, tprops['major'])    
    rotated, angles = rotate(scaled, tprops['angle'])
    aligned = crop_and_align(rotated, angles)
    result = np.mean(aligned , axis=0)
    return result, signaltonoise(result)



