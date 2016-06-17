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
import numpy.ma as ma
import pylab
import pyfits
import glob
import os
from PIL import Image
import pyfits
import numpy as np
import pylab as py
import img_scale


def normalizeAux(x,x1,x2,y1,y2):

    y = ( (x-x1)*(y2-y1)/(x2-x1) ) + y1
    return y

def normalize(data):
	maxv = data.max()
	minv = data.min()
	height, width = data.shape
	for i in xrange(0,height):
		for j in xrange(0,width):
			data[i][j] = normalizeAux(data[i][j], 100.0, 0.0, maxv, minv)

	return data

def cropAux (originalData):
	shape = originalData.shape
	if len(shape) == 4:
		originalData = originalData[0,0]

	height = len(originalData)
	width = len(originalData[0])
	data = np.ndarray(shape=(height,width), dtype=float)

	prom = originalData.mean(axis=1, dtype=float).mean()
	total = 0
	ntotal = 0

	minFilter = prom
	
	plusFilter = ma.masked_array(originalData, mask= originalData <= minFilter)
	total+= np.sum(plusFilter.compressed())
	ntotal+= plusFilter.compressed().shape[0]
	data = plusFilter.filled(0)

	#for i in xrange(0,height):
	#	for j in xrange(0,width):
	#		if originalData[i,j] > minFilter:
	#			data[i,j] = originalData[i,j]
	#			total+=originalData[i,j]
	#			ntotal+=1
	#		else:
	#			data[i,j] = 0

	minFilter = total/ntotal

	maxv = -999999999999
	# se suman al pixel actual los valores de todos los pixeles contiguos y se guarda el maximo valor
	

	for i in xrange(0,height):
		for j in xrange(0,width):
			if data[i,j] <= minFilter:
				continue
			tmp = minFilter
			if i > 0 and data[i-1,j] > minFilter: #abajo
				tmp += data[i-1,j]
			if i < height-1 and data[i+1,j] > minFilter: #arriba
				tmp += data[i+1,j]
			if j > 0 and data[i,j-1] > minFilter: #izquierda
				tmp += data[i,j-1]
			if j < width-1 and data[i,j+1] > minFilter: #derecha
				tmp += data[i,j+1]
			if i > 0 and j > 0 and data[i-1,j-1] > minFilter: #abajo-izquierda
				tmp += data[i-1,j-1]
			if i < height-1 and j > 0 and data[i+1,j-1] > minFilter: #arriba-izquierda
				tmp += data[i+1,j-1]
			if i > 0 and j < width-1 and data[i-1,j+1] > minFilter: #abajo-derecha
				tmp += data[i-1,j+1]
			if i < height-1 and j < width-1 and data[i+1,j+1] > minFilter: #arriba-derecha
				tmp += data[i+1,j+1]
			data[i,j] = tmp #+= tmp
			if data[i,j] > maxv:
				maxv = data[i,j]

	maxy = maxx = -1
	miny = minx = width+height

	total = 0
	ntotal = 0
	
	# al grupo que contiene el pixel con el mayor valor, se le setea el mayor valor
	for i in xrange(height-1,-1,-1):
		for j in xrange(width-1,-1,-1):
			if data[i,j] < maxv:
				continue
			
			total+=originalData[i,j]
			ntotal+=1

			if i > 0 and data[i-1,j] > minFilter: #abajo
				data[i-1,j] = maxv
			if i < height-1 and data[i+1,j] > minFilter: #arriba
				data[i+1,j] = maxv
			if j > 0 and data[i,j-1] > minFilter: #izquierda
				data[i,j-1] = maxv
			if j < width-1 and data[i,j+1] > minFilter: #derecha
				data[i,j+1] = maxv
			if i > 0 and j > 0 and data[i-1,j-1] > minFilter: #abajo-izquierda
				data[i-1,j-1] = maxv
			if i < height-1 and j > 0 and data[i+1,j-1] > minFilter: #arriba-izquierda
				data[i+1,j-1] = maxv
			if i > 0 and j < width-1 and data[i-1,j+1] > minFilter: #abajo-derecha
				data[i-1,j+1] = maxv
			if i < height-1 and j < width-1 and data[i+1,j+1] > minFilter: #arriba-derecha
				data[i+1,j+1] = maxv
			if i < miny:
				miny = i
			if i > maxy:
				maxy = i
			if j < minx:
				minx = j
			if j > maxx:
				maxx = j

	minFilter = total/ntotal

	maxy = maxx = -1
	miny = minx = width+height

	for i in xrange(0,height):
		for j in xrange(0,width):
			if data[i,j] < maxv or originalData[i,j] < minFilter:
				data[i, j] = 0
			else:
				data[i, j] = originalData[i,j]

				if i < miny:
					miny = i
				if i > maxy:
					maxy = i
				if j < minx:
					minx = j
				if j > maxx:
					maxx = j

	
	originalData = originalData[miny:maxy, minx:maxx]
	newHeight, newWidth = originalData.shape

	for i in xrange(0,newHeight):
		for j in xrange(0,newWidth):
			if originalData[i,j] > 0 and i == 0 and j == 0 and originalData[i+1,j] == 0 and originalData[i,j+1] == 0 and originalData[i+1,j+1] == 0:
				originalData[i,j] = 0
			elif originalData[i,j] > 0 and i == (newHeight-1) and j == (newWidth-1) and originalData[i-1,j] == 0 and originalData[i,j-1] == 0 and originalData[i-1,j-1] == 0:
				originalData[i,j] = 0
			elif originalData[i,j] > 0 and i == 0 and j == (newWidth-1) and originalData[i+1,j] == 0 and originalData[i,j-1] == 0 and originalData[i+1,j-1] == 0:
				originalData[i,j] = 0
			elif originalData[i,j] > 0 and i == (newHeight-1) and j == 0 and originalData[i-1,j] == 0 and originalData[i,j+1] == 0 and originalData[i-1,j+1] == 0:
				originalData[i,j] = 0
			elif originalData[i,j] > 0 and i == (newHeight-1) and j > 0 and j < (newWidth-1) and originalData[i-1,j] == 0 and originalData[i,j+1] == 0 and originalData[i,j-1] == 0 and originalData[i-1,j+1] == 0 and originalData[i-1,j-1] == 0:
				originalData[i,j] = 0
			elif originalData[i,j] > 0 and i < (newHeight-1) and i > 0 and j == 0 and originalData[i-1,j] == 0 and originalData[i-1,j+1] == 0 and originalData[i,j+1] == 0 and originalData[i+1,j+1] == 0 and originalData[i+1,j] == 0:
				originalData[i,j] = 0
			elif originalData[i,j] > 0 and i == 0 and j > 0 and j < (newWidth-1) and originalData[i,j-1] == 0 and originalData[i,j+1] == 0 and originalData[i+1,j+1] == 0 and originalData[i+1,j] == 0 and originalData[i+1,j-1] == 0:
				originalData[i,j] = 0
			elif originalData[i,j] > 0 and i > 0 and i < (newHeight-1) and j == (newWidth-1) and originalData[i+1,j] == 0 and originalData[i+1,j-1] == 0 and originalData[i,j-1] == 0 and originalData[i-1,j-1] == 0 and originalData[i-1,j] == 0:
				originalData[i,j] = 0
			elif originalData[i,j] > 0 and i > 0 and i < (newHeight-1) and j < (newWidth-1) and j > 0 and originalData[i+1,j] == 0 and originalData[i+1,j-1] == 0 and originalData[i,j-1] == 0 and originalData[i-1,j-1] == 0 and originalData[i-1,j] == 0 and originalData[i-1,j+1] == 0 and originalData[i,j+1] == 0 and originalData[i+1,j+1] == 0:
				originalData[i,j] = 0
			elif originalData[i,j] > 0:
				if i < miny:
					miny = i
				if i > maxy:
					maxy = i
				if j < minx:
					minx = j
				if j > maxx:
					maxx = j
	
	newData = originalData[miny:maxy, minx:maxx]
	newHeight, newWidth = newData.shape

	border = []
	for i in xrange(0,newHeight): 
		line = []
		for j in xrange(0,newWidth):
			if newData[i,j] > 0:
				line.append((i,j))

		if len(line) > 0:
			border.append(line[0])
		if len(line) > 1:
			border.append(line[-1])

	newData = normalize(newData)

	return border, newData

def minValue(data):

	x,y = data.shape
	minvalue = -99999999 

	mindata = data.min()

	if mindata < minvalue:
		minvalue = mindata
	
	#for i in xrange(0,x):
	#	for j in xrange(0,y):
	#		if data[i][j] < minvalue:
	#			minvalue = data[i][j]

	return minvalue

def manualCrop(data,E1,E2,E3,E4):

	newdata = data[int(E3):int(E4), int(E1):int(E2)]
	return newdata


def crop(inputDir, outputDir):
	data = glob.glob(inputDir+'/*.fits')
	borders = []
	# - - - - - -  CREA EL SUBDIRECTORIO PARA GUARDAR LAS FOTOS PNG - - - - - - 
	dir_png = outputDir+'/PNG_Images' 
	if not os.path.isdir(dir_png):
		os.makedirs(dir_png)
	

	for i in xrange(0,len(data)):

		name = data[i].split('/')[-1]#.split('.')[0]
		image = fits.open(data[i], ignore_missing_end = True )[0].data
		if isinstance(image, list):
			image = image[0]

		print "Crop: "+'/Img_0_'+str(i)+'.fits'

		fits.writeto(outputDir+'/Img_0_'+str(i)+'.fits',image, clobber=True)	
		border, image = cropAux(image)
		fits.writeto(outputDir+'/Img_1_'+str(i)+'.fits',image, clobber=True)
		borders.append(border)
		# - - - - - -  CREA Y GUARDA LAS FOTOS PNG EN EL SUBDIRECTORIO - - - - - - 
		# j_img = pyfits.getdata(outputDir+'/Img_1_'+str(i)+'.fits')
		# img = np.zeros((j_img.shape[0], j_img.shape[1]), dtype=float)
		# img[:,:] = img_scale.sqrt(j_img, scale_min=0, scale_max=10000)
		# py.clf()
		# py.imshow(img, aspect='equal')
		# py.title('Crop Img_1_'+str(i))
		# py.savefig(dir_png+'/Img_1_'+str(i)+'.png')
		# img = Image.open(dir_png+'/Img_1_'+str(i)+'.png')
		# img.show()

	print "Done."
	return borders

def cropManual(inputDir, outputDir,E1,E2,E3,E4):
	data = glob.glob(inputDir+'/*.fits')
	borders = []
	# - - - - - -  CREA EL SUBDIRECTORIO PARA GUARDAR LAS FOTOS PNG - - - - - - 
	dir_png = outputDir+'/PNG_Images' 
	if not os.path.isdir(dir_png):
		os.makedirs(dir_png)
	for i in xrange(0,len(data)):

		name = data[i].split('/')[-1].split('.')[0]
		image = fits.getdata(data[i])
		if isinstance(image, list):
			image = image[0]

		print "Crop: "+'/Img_0_'+str(i)+'.fits'

		fits.writeto(outputDir+'/Img_0_'+str(i)+'.fits',image, clobber=True)	
		image = manualCrop(image,E1,E2,E3,E4)
		fits.writeto(outputDir+'/Img_1_'+str(i)+'.fits',image, clobber=True)

	print "Done."


