import sys
sys.path.append('../../')

from mayavi import mlab

import acalib.vo.workspace as ws
import time
from scipy import ndimage 
import numpy as np
import matplotlib.pyplot as plt
import timeit
import cProfile
import acalib.clumps.bubbleClumps as cl
import acalib.clumps.gaussClumps as gl
import matplotlib
matplotlib.use('WxAgg')
matplotlib.interactive(True)
import matplotlib.pyplot as plt


binpath='../../bindata/fits/cubes/'
#ws.import_file(binpath+"M100line.image.fits")
ws.import_file(binpath+"Orion.methanol.cbc.contsub.image.fits")
#ws.import_file(binpath+"Antennae_North.CO3_2Line.Clean.pcal1.image.fits")
#ws.import_file(binpath+"Antennae_South.CO3_2Line.Clean.pcal1.image.fits")
#ws.import_file(binpath+"CenA.CO2_1Line.Clean.image.fits")

elm=ws.elements()
#cube=elm['M100line.image-0']
cube=elm['Orion.methanol.cbc.contsub.image-0']
#cube=elm['Antennae_North.CO3_2Line.Clean.pcal1.image-0']
#cube=elm['Antennae_South.CO3_2Line.Clean.pcal1.image-0']
#cube=elm['CenA.CO2_1Line.Clean.image-0']
spar=cube.standarize()
bc=cl.BubbleClumps()
# use_meta not implemented yet, so compute parameters to use
pixbsize=cube.meta['BMIN']/abs(cube.meta['CDELT1'])
print "beam size in pixels =",pixbsize
bc.par['FWHMBEAM']=pixbsize

gc=gl.GaussClumps()
gc.par['FWHMBEAM']=pixbsize

# Fitme :S
clist=gc.fit(cube,verbose=True)

porc=[]
elem=[]
snr=[]
tim=[]


telem=float(cube.count())
maxporc=0.01
samples=maxporc*telem
print "Maximum Bubbles",int(samples),"= ",maxporc*100,"%" 
bc.par['MAXBUB']=int(samples)
snrlimit=0.5
print "SNR = ",snrlimit
bc.par['SNRLIMIT']=snrlimit
bc.fit(cube,verbose=True)


# Fitering and Thresholding
rms=cube.estimate_rms()
zcube=np.nan_to_num(cube.data)
newcube=cube.empty_like()
newcube.data=ndimage.gaussian_filter(zcube,[bc.ss,bc.sb,bc.sb])
newcube.data[newcube.data<=rms +  snrlimit*rms]=0

#cube.stacked_show()
#newcube.stacked_show()
#bc.syn.stacked_show()

#cube.volume_show()
#newcube.volume_show()
#bc.syn.volume_show()

#cube.contour_show()
newcube.contour_show()
gc.syn.contour_show()
bc.syn.contour_show()


emptycube=cube.empty_like()
plt.subplot(2,4,1)
plt.imshow(cube.stack())
plt.subplot(2,4,2)
plt.imshow(newcube.stack())
plt.subplot(2,4,3)
plt.imshow(gc.syn.stack())
plt.subplot(2,4,4)
plt.imshow(bc.syn.stack())
plt.subplot(2,4,5)
plt.imshow(emptycube.stack())
plt.subplot(2,4,6)
plt.imshow(cube.stack() - newcube.stack())
plt.subplot(2,4,7)
plt.imshow(gc.data.stack())
plt.subplot(2,4,8)
plt.imshow(bc.residual.stack())
plt.show()
plt.pause(100)
