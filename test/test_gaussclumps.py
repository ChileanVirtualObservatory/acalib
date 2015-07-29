import core.workspace as ws
import numpy as np
import matplotlib.pyplot as plt
import timeit
import cProfile
import clumps.gaussClumps as gclumps

ws.import_file("fits/M100line.image.fits")

elm=ws.elements()
cube=elm['M100line.image-0']
gc=gclumps.GaussClumps()
# use_meta not implemented yet, so compute parameters to use
pixbsize=cube.meta['BMIN']/abs(cube.meta['CDELT1'])
print "beam size in pixels =",pixbsize
gc.par['FWHMBEAM']=pixbsize

# Fitme :S
gc.fit(cube,verbose=True)





