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
bc=cl.BubbleClumps()
# use_meta not implemented yet, so compute parameters to use
pixbsize=cube.meta['BMIN']/abs(cube.meta['CDELT1'])
print "beam size in pixels =",pixbsize
bc.par['FWHMBEAM']=pixbsize

total=cube.flux()
telem=float(cube.count())
maxporc=0.01
samples=maxporc*telem
print "Maximum Bubbles",int(samples),"= ",maxporc*100,"%" 
bc.par['MAXBUB']=int(samples)
snrlimit=0.5
print "SNR = ",snrlimit
bc.par['SNRLIMIT']=snrlimit
bc.fit(cube,verbose=True)

gp.volume(bc.syn)
gp.contour(bc.syn)

clust=bc.clustering(10.0,method='dbscan')
fig = plt.figure("DBSCAN")
bc.draw_cluster(fig,clust)

clust=bc.clustering(0.8,method='affinity_propagation')
fig = plt.figure("AFFINITY PROPAGATION")
bc.draw_cluster(fig,clust)

clust=bc.clustering(5,method='kmeans')
fig = plt.figure("KMEANS")
bc.draw_cluster(fig,clust)

clust=bc.clustering(5,method='spectral')
fig = plt.figure("SPECTRAL")
bc.draw_cluster(fig,clust)

plt.show()


