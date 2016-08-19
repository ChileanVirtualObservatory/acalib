

import props as pr
import transforms as tr

import matplotlib.pyplot as plt
import numpy as np


def stacking(template_path, data_path):
	tprops = pr.fits_props(template_path)

	scaled = tr.scale(data_path, tprops['major'])

	rotated = tr.rotate(scaled, tprops['angle'])

	aligned = tr.cropAndAlign(rotated, tprops['angle'])

	result = np.mean(aligned, axis = 0)

	plt.imshow(result)
	plt.show

