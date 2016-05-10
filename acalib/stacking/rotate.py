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
import glob
import math
from PIL import Image
import pyfits
import numpy as np
import pylab as py
import img_scale

def rotate_coords(x, y, theta, ox, oy):
	s, c = np.sin(theta), np.cos(theta)
	x, y = np.asarray(x) - ox, np.asarray(y) - oy
	return x * c - y * s + ox, x * s + y * c + oy

def rotate_image(src, theta, fill=0):
	# theta = -theta*np.pi/180
	oy = len(src)
	ox = len(src[0])

	sh, sw = src.shape

	cx, cy = rotate_coords([0, sw, sw, 0], [0, 0, sh, sh], theta, ox, oy)

	dw, dh = (int(np.ceil(c.max() - c.min())) for c in (cx, cy))

	dx, dy = np.meshgrid(np.arange(dw), np.arange(dh))

	sx, sy = rotate_coords(dx + cx.min(), dy + cy.min(), -theta, ox, oy)

	sx, sy = sx.round().astype(int), sy.round().astype(int)

	mask = (0 <= sx) & (sx < sw) & (0 <= sy) & (sy < sh)

	dest = np.empty(shape=(dh, dw), dtype=src.dtype)

	dest[dy[mask], dx[mask]] = src[sy[mask], sx[mask]]

	dest[dy[~mask], dx[~mask]] = fill

	return dest

def fartestPoints(border):
	dist = []
	for i in xrange(0,len(border)):
		for j in xrange(0,len(border)):
			y1,x1 = border[i]
			y2,x2 = border[j]
			d = (x1-x2)**2 + (y1-y2)**2
			d = math.sqrt(d)
			dist.append((d,border[i],border[j]))
	maxd = 0
	points = 0
	for i in xrange(0,len(dist)):
		if dist[i][0] > maxd:
			points = (dist[i][1],dist[i][2])
			maxd = dist[i][0]

	return points

def theta(points):
	print points
	y1,x1 = points[0]
	y2,x2 = points[1]
	theta = (y2-y1+.0)/(x2-x1+.0)
	theta = math.atan(theta+(math.pi*9)/2)
	return theta

def rotate(outputDir, border):
	maxHeight = 0
	maxWidth = 0
	data = sorted(glob.glob(outputDir+'/Img_1_*.fits'))
	dir_png = outputDir+'/PNG_Images' 
	for i in xrange(0,len(data)):
		image = fits.getdata(data[i])
		points = fartestPoints(border[i])
		image = rotate_image(image,theta(points))

		print "Rotate: "+'/Img_1_'+str(i)+'.fits',

		fits.writeto(outputDir+'/Img_2_'+str(i)+'.fits',image, clobber=True)

		# j_img = pyfits.getdata(outputDir+'/Img_2_'+str(i)+'.fits')
		# img = np.zeros((j_img.shape[0], j_img.shape[1]), dtype=float)
		# img[:,:] = img_scale.sqrt(j_img, scale_min=0, scale_max=10000)
		# py.clf()
		# py.imshow(img, aspect='equal')
		# py.title('Rotate Img_2_'+str(i))
		# py.savefig(dir_png+'/Img_2_'+str(i)+'.png')
		# img = Image.open(dir_png+'/Img_2_'+str(i)+'.png')
		# img.show()
		
		print "Done."

		h,w = image.shape
		if h > maxHeight:
			maxHeight = h
		if w > maxWidth:
			maxWidth = w

	return (maxHeight,maxWidth)

def rotateManual(outputDir):

	data = sorted(glob.glob(outputDir+'/Img_1_*.fits'))
	# dir_png = outputDir+'/PNG_Images' 
	for i in xrange(0,len(data)):
		image = fits.getdata(data[i])
		image = rotate_image(image,math.atan(theta+(math.pi*9)/2))

		print "Rotate: "+'/Img_1_'+str(i)+'.fits',

		fits.writeto(outputDir+'/Img_2_'+str(i)+'.fits',image, clobber=True)

		# j_img = pyfits.getdata(outputDir+'/Img_2_'+str(i)+'.fits')
		# img = np.zeros((j_img.shape[0], j_img.shape[1]), dtype=float)
		# img[:,:] = img_scale.sqrt(j_img, scale_min=0, scale_max=10000)
		# py.clf()
		# py.imshow(img, aspect='equal')
		# py.title('Rotate Img_2_'+str(i))
		# py.savefig(dir_png+'/Img_2_'+str(i)+'.png')
		# img = Image.open(dir_png+'/Img_2_'+str(i)+'.png')
		# img.show()
		
		print "Done."
