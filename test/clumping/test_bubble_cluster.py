import sys
sys.path.append('../../')


import acalib.vo.workspace as ws
import time
import numpy as np
import matplotlib.pyplot as plt

import timeit
import cProfile
import acalib.clumps.bubbleClumps as cl
import matplotlib.pyplot as plt



binpath='../../bindata/fits/cubes/'
#ws.import_file(binpath+"M100line.image.fits")
ws.import_file(binpath+"Orion.methanol.cbc.contsub.image.fits")
#ws.import_file(binpath+"Boom.cm.cln.fits")
#ws.import_file(binpath+"Antennae_North.CO3_2Line.Clean.pcal1.image.fits")
#ws.import_file(binpath+"Antennae_South.CO3_2Line.Clean.pcal1.image.fits")

elm=ws.elements()
#cube=elm['M100line.image-0']
cube=elm['Orion.methanol.cbc.contsub.image-0']
#cube=elm['Boom.cm.cln-0']
#cube=elm['Antennae_North.CO3_2Line.Clean.pcal1.image-0']
#cube=elm['Antennae_South.CO3_2Line.Clean.pcal1.image-0']
spar=cube.standarize()

spar=cube.standarize()
gc=cl.BubbleClumps()
# use_meta not implemented yet, so compute parameters to use
pixbsize=cube.meta['BMIN']/abs(cube.meta['CDELT1'])
print "beam size in pixels =",pixbsize
gc.par['FWHMBEAM']=pixbsize

total=cube.flux()
telem=float(cube.count())
maxporc=0.01
samples=maxporc*telem
print "Maximum Bubbles",int(samples),"= ",maxporc*100,"%" 
gc.par['MAXBUB']=int(samples)
snrlimit=0.5
print "SNR = ",snrlimit
gc.par['SNRLIMIT']=snrlimit
gc.fit(cube,verbose=True)

cube.volume_show()
cube.contour_show()
gc.syn.volume_show()
gc.syn.contour_show()

#gc.test_clustering()
gc.selected_clusters(5)
#gc.reasonable_cluster()

