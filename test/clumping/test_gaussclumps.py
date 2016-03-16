import sys
sys.path.append('../../')

from acalib import acontainer as ac
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

gc.syn.stacked_show()
gc.syn.volume_show()
gc.syn.contour_show()

gc.data.stacked_show()
gc.data.volume_show()
gc.data.contour_show()


