import sys
sys.path.append('../../')

from acalib import acontainer as ac
from acalib.io import graph as gp
import numpy as np
import matplotlib.pyplot as plt
import timeit
import cProfile
import acalib.clumps.gaussClumps as gclumps
import matplotlib.pyplot as plt

cont=ac.AContainer()
cont.load(sys.argv[1])
cube=cont.primary
spar=cube.standarize()
gc=gclumps.GaussClumps()
# use_meta not implemented yet, so compute parameters to use
pixbsize=cube.meta['BMIN']/abs(cube.meta['CDELT1'])
print "beam size in pixels =",pixbsize
gc.par['FWHMBEAM']=pixbsize

# Fitme :S
clist=gc.fit(cube,verbose=True)
print clist

gp.velocity(cube)
gp.stacked(cube)
gp.volume(cube)
gp.contour(cube)

gp.stacked(gc.syn)
gp.volume(gc.syn)
gp.contour(gc.syn)

gp.stacked(gc.syn)
gp.volume(gc.syn)
gp.contour(gc.syn)



