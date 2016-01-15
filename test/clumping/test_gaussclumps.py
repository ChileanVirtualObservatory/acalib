import sys
sys.path.append('../../')

import acalib.vo.workspace as ws
import numpy as np
import matplotlib.pyplot as plt
import timeit
import cProfile
import acalib.clumps.gaussClumps as gclumps
import matplotlib.pyplot as plt

binpath='../../bindata/fits/cubes/'
ws.import_file(binpath+"M100line.image.fits")
#ws.import_file(binpath+"Orion.methanol.cbc.contsub.image.fits")

elm=ws.elements()
cube=elm['M100line.image-0']
#cube=elm['Orion.methanol.cbc.contsub.image-0']
spar=cube.standarize()
gc=gclumps.GaussClumps()
# use_meta not implemented yet, so compute parameters to use
pixbsize=cube.meta['BMIN']/abs(cube.meta['CDELT1'])
print "beam size in pixels =",pixbsize
gc.par['FWHMBEAM']=pixbsize

# Fitme :S
clist=gc.fit(cube,verbose=True)
for clump in clist:
   print("a=",clump[0],"b=",clump[1],"pos",(clump[2],clump[4],clump[7]),"std",(clump[3],clump[5],clump[8]),"ang",clump[6])

#cube.stacked_show()
#cube.volume_show()
#cube.contour_show()

gc.syn.stacked_show()
gc.syn.volume_show()
gc.syn.contour_show()


gc.data.stacked_show()
gc.data.volume_show()
gc.data.contour_show()

