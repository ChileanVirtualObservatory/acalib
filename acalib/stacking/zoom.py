import numpy as np
from scipy.ndimage.interpolation import zoom
from astropy.io import fits
import os
import glob
import matplotlib.pyplot as plt

def magnitude(inputDir, outputDir, x):
	data = glob.glob(inputDir+'/*.fits')


	for i in xrange(0,len(data)):
		name = data[i].split('/')[1]
		
		image = fits.open(data[i], ignore_missing_end = True )[0].data
		new_image = zoom(image, x)
		fits.writeto(outputDir+'/'+name,new_image, clobber=True)


def images(inputDir, outputDir):
	data = glob.glob(inputDir+'/*.fits')
	for i in xrange(0,len(data)):
		image = fits.open(data[i], ignore_missing_end = True )[0].data
		name = data[i].split('/')[1].split('.')[0]
		plt.imshow(image)
		plt.savefig(outputDir+'/'+name)





