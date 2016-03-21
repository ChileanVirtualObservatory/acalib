import sys
sys.path.append('../../')

import time
import numpy as np
import matplotlib.pyplot as plt

import timeit
import cProfile
import acalib.clumps.bubbleClumps as cl
from acalib import acontainer as ac
from acalib.io import graph as gp

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

gp.volume(cube)
gp.countour(cube)
gp.volume(gc.syn)
gp.countour(gc.syn)

#gc.test_clustering()
gc.selected_clusters(5)
#gc.reasonable_cluster()

