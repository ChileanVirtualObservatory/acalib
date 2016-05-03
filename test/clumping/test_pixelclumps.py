import sys
sys.path.append('../../')

from acalib import acontainer as ac
from acalib.io import graph as gp
import numpy as np
import matplotlib.pyplot as plt
import timeit
import cProfile
import acalib.clumps.pixelClumps as pclumps
import matplotlib.pyplot as plt

cont=ac.AContainer()
cont.load(sys.argv[1])
cube=cont.primary
spar=cube.standarize()
pc=pclumps.PixelClumps()
# use_meta not implemented yet, so compute parameters to use
pixbsize=cube.meta['BMIN']/abs(cube.meta['CDELT1'])
print "beam size in pixels =",pixbsize
pc.par['FWHMBEAM']=pixbsize
pc.par['SNRLIMIT']=1.0

# Fitme :S
clist=pc.fit(cube,verbose=True)
print clist

#gp.velocity(cube)
#gp.stacked(cube)
#gp.volume(cube)
#gp.contour(cube)

#gp.stacked(gc.syn)
#gp.volume(gc.syn)
#gp.contour(gc.syn)

#gp.stacked(gc.syn)
#gp.volume(gc.syn)
#gp.contour(gc.syn)



