#This file is part of ChiVO, the Chilean Virtual Observatory
#A project sponsored by FONDEF (D11I1060)
#Copyright (C) 2015 Universidad Tecnica Federico Santa Maria Mauricio Solar
#                                                            Marcelo Mendoza
#                   Universidad de Chile                     Diego Mardones
#                   Pontificia Universidad Catolica          Karim Pichara
#                   Universidad de Concepcion                Ricardo Contreras
#                   Universidad de Santiago                  Victor Parada
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

from astropy.io import fits
import numpy as np
import math
import glob
from PIL import Image
import pyfits
import numpy as np
import pylab as py
import img_scale

import matplotlib as plt

def stacking(outputDir,maxSize):

	data = sorted(glob.glob(outputDir+'/Img_4_*.fits'))
	dir_png = outputDir+'/PNG_Images'	

	newdata = np.zeros((maxSize[0],maxSize[1]))
	datasum = np.zeros((maxSize[0],maxSize[1]))

	for x in xrange(0,len(data)):
		img = fits.getdata(data[x])
		for i in xrange(0,len(img)):
			for j in xrange(0,len(img[i])):
				newdata[i][j] += img[i][j]
				if newdata[i][j] > 0:
					datasum[i][j] += 1 

	# for i in xrange(0,len(newdata)):
	# 	for j in xrange(0,len(newdata[0])):
	# 			if datasum[i][j] != 0:
	# 				newdata[i][j] /= datasum[i][j]

	

	fits.writeto(outputDir+'/Img_5.fits', newdata, clobber=True)
	# j_img = pyfits.getdata(outputDir+'/Img_5.fits')
	# img = np.zeros((j_img.shape[0], j_img.shape[1]), dtype=float)
	# img[:,:] = img_scale.sqrt(j_img, scale_min=0, scale_max=10000)
	# py.clf()
	# py.imshow(img, aspect='equal')
	# py.title('Stack Img_5')
	# py.savefig(dir_png+'/Img_5.png')
	# img = Image.open(dir_png+'/Img_5.png')
	# img.show()


	print "Stack: Done."	

def stackingManual(outputDir):

	data = sorted(glob.glob(outputDir+'/Img_1_*.fits'))
	# dir_png = outputDir+'/PNG_Images'	
	img = fits.getdata(data[0])
	h,w = img.shape
	newdata = np.zeros((h,w))
	datasum = np.zeros((h,w))

	for x in xrange(0,len(data)):
		img = fits.getdata(data[x])
		for i in xrange(0,len(img)):
			for j in xrange(0,len(img[i])):
				newdata[i][j] += img[i][j]
				if newdata[i][j] > 0:
					datasum[i][j] += 1 

	# for i in xrange(0,len(newdata)):
	# 	for j in xrange(0,len(newdata[0])):
	# 			if datasum[i][j] != 0:
	# 				newdata[i][j] /= datasum[i][j]
	

	fits.writeto(outputDir+'/Img_2.fits', newdata, clobber=True)
	# j_img = pyfits.getdata(outputDir+'/Img_5.fits')
	# img = np.zeros((j_img.shape[0], j_img.shape[1]), dtype=float)
	# img[:,:] = img_scale.sqrt(j_img, scale_min=0, scale_max=10000)
	# py.clf()
	# py.imshow(img, aspect='equal')
	# py.title('Stack Img_5')
	# py.savefig(dir_png+'/Img_5.png')
	# img = Image.open(dir_png+'/Img_5.png')
	# img.show()


	print "Stack: Done."	
