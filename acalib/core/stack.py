

import props as pr
import transforms as tr

import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import signaltonoise

from utils import moment2


def stacking(template_path, data_path):
	tprops = pr.fits_props(template_path)

	scaled = tr.scale(data_path, tprops['major'])

	rotated, angles = tr.rotate(scaled, tprops['angle'])

	aligned = tr.cropAndAlign(rotated, angles)

	result = np.mean(aligned, axis = 0)
	return result, signaltonoise(result), moment2(result) 	
	#plt.imshow(result)
	#plt.show()



