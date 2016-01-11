import acalib.core.workspace as ws
import time
import numpy as np
import matplotlib.pyplot as plt
import timeit
import cProfile
import acalib.clumps.bubbleClumps as cl
import matplotlib.pyplot as plt

ws.import_file("fits/M100line.image.fits")
#ws.import_file("fits/Orion.methanol.cbc.contsub.image.fits")
#ws.import_file("fits/Antennae_North.CO3_2Line.Clean.pcal1.image.fits")
#ws.import_file("fits/Antennae_South.CO3_2Line.Clean.pcal1.image.fits")
#ws.import_file("fits/CenA.CO2_1Line.Clean.image.fits")

elm=ws.elements()
cube=elm['M100line.image-0']
#cube=elm['Orion.methanol.cbc.contsub.image-0']
#cube=elm['Antennae_North.CO3_2Line.Clean.pcal1.image-0']
#cube=elm['Antennae_South.CO3_2Line.Clean.pcal1.image-0']
#cube=elm['CenA.CO2_1Line.Clean.image-0']
spar=cube.standarize()
gc=cl.BubbleClumps()
# use_meta not implemented yet, so compute parameters to use
pixbsize=cube.meta['BMIN']/abs(cube.meta['CDELT1'])
print "beam size in pixels =",pixbsize
gc.par['FWHMBEAM']=pixbsize

total=cube.sum()
telem=float(cube.count())
maxporc=0.001
samples=maxporc*telem
print "Maximum Bubbles",int(samples),"= ",maxporc*100,"%" 
gc.par['MAXBUB']=int(samples)
snrlimit=0.5
print "SNR = ",snrlimit
gc.par['SNRLIMIT']=snrlimit
gc.fit(cube,verbose=True)
#gc.test_clustering()
gc.cluster_candidates(method='affinity_propagation')

