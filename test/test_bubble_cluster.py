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
#ws.import_file("fits/calibrated.ms.line.spw0.source15.image.fits")
#ws.import_file("fits/HOT2_EI__e_03_TE_source.6.13CO.image.pbcor.fits")
#ws.import_file("fits/Boom.cm.cln.fits")


elm=ws.elements()
cube=elm['M100line.image-0']
#cube=elm['Orion.methanol.cbc.contsub.image-0']
#cube=elm['Antennae_North.CO3_2Line.Clean.pcal1.image-0']
#cube=elm['Antennae_South.CO3_2Line.Clean.pcal1.image-0']
#cube=elm['calibrated.ms.line.spw0.source15.image-0']
#cube=elm['HOT2_EI__e_03_TE_source.6.13CO.image.pbcor-0']
#cube=elm['Boom.cm.cln-0']

spar=cube.standarize()
gc=cl.BubbleClumps()
# use_meta not implemented yet, so compute parameters to use
pixbsize=cube.meta['BMIN']/abs(cube.meta['CDELT1'])
print "beam size in pixels =",pixbsize
gc.par['FWHMBEAM']=pixbsize

total=cube.sum()
telem=float(cube.count())
maxporc=0.01
samples=maxporc*telem
print "Maximum Bubbles",int(samples),"= ",maxporc*100,"%" 
gc.par['MAXBUB']=int(samples)
snrlimit=0.5
print "SNR = ",snrlimit
gc.par['SNRLIMIT']=snrlimit
gc.fit(cube,verbose=True)
#gc.test_clustering()
gc.selected_clusters(5)
#gc.reasonable_cluster()

