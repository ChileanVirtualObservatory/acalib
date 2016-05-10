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
import os

import crop as cr
import rotate as rt
import scale as sc
import align as al
import stacking as st
import compare as cm


def stack(inputDir,outputDir,E1,E2,E3,E4):
	print 'Input Directory:', inputDir
	print 'Output Directory:', outputDir
	for the_file in os.listdir(outputDir):
		file_path = os.path.join(outputDir, the_file)
		try:
			if os.path.isfile(file_path):
				os.unlink(file_path)
		except Exception, e:
			print e
	if (E1 or E2 or E3 or E4) == None:
		print 'A'
		border = cr.crop(inputDir, outputDir)
		maxSize = rt.rotate(outputDir, border)
		maxSize = sc.scale(outputDir, maxSize)
		al.align(outputDir, maxSize)
		st.stacking(outputDir, maxSize)
	elif(E1 and E2 and E3 and E4) != None:
		print 'B'
		cr.cropManual(inputDir,outputDir,E1,E2,E3,E4)
		st.stackingManual(outputDir)

def aux(inputDir,outputDir):
	data = glob.glob(inputDir+'/*.fits')
	print data
	print "1", data[0].split('/')[-1].split('.')[0]
	image1 = fits.getdata(data[0])
	print "2", data[1].split('/')[-1].split('.')[0]
	image2 = fits.getdata(data[1])
	cm.compare(image1, image2)
	
