import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import signaltonoise
from acalib.io.fits import load_fits_to_cont

def stacking(template_data, data_cont):
	tprops = pr.fits_props(template_data)

	scaled = tr.scale(data_cont, tprops['major'])

	rotated, angles = tr.rotate(scaled, tprops['angle'])

	aligned = tr.cropAndAlign(rotated, angles)

	result = np.mean(aligned, axis = 0)




	return result, signaltonoise(result), moment2(result) 	
	#plt.imshow(result)
	#plt.show()



