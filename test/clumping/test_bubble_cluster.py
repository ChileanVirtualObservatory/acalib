import sys
sys.path.append('../../')

import time
import numpy as np
import matplotlib.pyplot as plt

import timeit
import cProfile
import acalib.clumps.bubbleClumps as cl
from acalib import acontainer as ac

import matplotlib.pyplot as plt


cont=ac.AContainer()
cont.load(sys.argv[1])
cube=cont.primary


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

