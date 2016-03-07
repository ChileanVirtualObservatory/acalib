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
import matplotlib
matplotlib.use('WxAgg')
matplotlib.interactive(True)
import matplotlib.pyplot as plt
from matplotlib import cm

import acalib.clumps.bubbleClumps as bclumps
import acalib.clumps.gaussClumps as gclumps
import acalib.clumps.fellWalker as fwalker


binpath='../../bindata/fits/cubes/'
#ws.import_file(binpath+"M100line.image.fits")
#ws.import_file(binpath+"Orion.methanol.cbc.contsub.image.fits")
#ws.import_file(binpath+"Boom.cm.cln.fits")
#ws.import_file(binpath+"Antennae_North.CO3_2Line.Clean.pcal1.image.fits")
#ws.import_file(binpath+"Antennae_South.CO3_2Line.Clean.pcal1.image.fits")
ws.import_file(binpath+"calibrated.ms.contsub.bin4.line.fits")
#ws.import_file(binpath+"calibrated.ms.contsub.bin4.line.image.fits")

elm=ws.elements()
#cube=elm['M100line.image-0']
#cube=elm['Orion.methanol.cbc.contsub.image-0']
#cube=elm['Boom.cm.cln-0']
#cube=elm['Antennae_North.CO3_2Line.Clean.pcal1.image-0']
#cube=elm['Antennae_South.CO3_2Line.Clean.pcal1.image-0']
cube=elm['calibrated.ms.contsub.bin4.line-0']
#cube=elm['calibrated.ms.contsub.bin4.line.image-0']
spar=cube.standarize()

# use_meta not implemented yet, so compute parameters to use
pixbsize=cube.meta['BMIN']/abs(cube.meta['CDELT1'])
snrlimit=1.5
print "[TEST] beam size in pixels =",pixbsize

# Bubble Detection
def BCProc():
   # Bubble Clump
   bc=bclumps.BubbleClumps()
   bc.par['FWHMBEAM']=pixbsize
   telem=float(cube.count())
   maxporc=0.01
   samples=maxporc*telem
#   print "[TEST] Maximum Bubbles",int(samples),"= ",maxporc*100,"%" 
   bc.par['MAXBUB']=int(samples)
   print "[TEST] SNR limit =",snrlimit
   bc.par['SNRLIMIT']=snrlimit
   bc.par['BSIZE']=1.0
   bc.fit(cube,verbose=True)
   
   #bc.fit(cube)
   print "Bubbles"
   print "elem/total ",float(len(bc.amplitudes))/float(cube.data.count())
   return bc.syn

# GaussClump Detection
def GCProc():
   gc=gclumps.GaussClumps()
   gc.par['FWHMBEAM']=pixbsize
   gc.par['THRESH']=snrlimit
   clist=gc.fit(cube)
   clist=gc.fit(cube,verbose=True)
   print "GaussClump"
   print "elem/total ",float(len(clist))/float(cube.data.count())
   return gc.syn

def GFProc():
# Fitering and Thresholding
   FWHM_TO_SIGMA = 1. / (8 * np.log(2))**0.5
   rms=cube.estimate_rms()
   zcube=np.nan_to_num(cube.data)
   newcube=cube.empty_like()
   
   sb=pixbsize*FWHM_TO_SIGMA
   ss=2*FWHM_TO_SIGMA

   newcube.data=ndimage.gaussian_filter(zcube,[ss,sb,sb])
   temp=cube.data[newcube.data>rms +  snrlimit*rms]
   temp2=newcube.data[newcube.data>rms +  snrlimit*rms]
   fact=(temp/temp2).min()
   newcube.data[newcube.data<=rms +  snrlimit*rms]=0
   newcube.data=fact*np.ma.masked_array(newcube.data,mask=cube.data.mask)
   print "Gaussian Filter"
   print "elem/total ",float(newcube.data[newcube.data>rms +  snrlimit*rms].count())/float(cube.data.count())
   return newcube

def PTProc():
   rms=cube.estimate_rms()
   newcube=cube.copy()
   newcube.data[newcube.data<=rms +  snrlimit*rms]=0
   newcube.data[newcube.data>rms +  snrlimit*rms]=newcube.data[newcube.data>rms +  snrlimit*rms] - rms
   print "PointThreshold"
   print "elem/total ",float(newcube.data[newcube.data>rms +  snrlimit*rms].count())/float(cube.data.count())
   return newcube


#def FWProc():
#   rms=cube.estimate_rms()
#   newcube=cube.copy()
#   newcube.data[newcube.data<=rms +  snrlimit*rms]=0
#   newcube.data[newcube.data>rms +  snrlimit*rms]=newcube.data[newcube.data>rms +  snrlimit*rms] - rms
#   print "PointThreshold"
#   print "elem/total ",float(newcube.data[newcube.data>rms +  snrlimit*rms].count())/float(cube.data.count())
#   return newcube

# Show Original Cube

cmap=cm.jet
cubelist=[GCProc(),GFProc(),PTProc(), BCProc()]
names=["GaussClumps","GaussianFilter","PointThreshold","BubbleClumps"]
n=len(cubelist)
# Compute residuals
res=[]
max_re=0
for i in range(0,n):
   cb=cubelist[i]
   re=cube.stack() - cb.stack()
   remax=re.max() 
   if remax > max_re:
      max_re=remax
   res.append(re)
   #cb.volume_show()
   #cb.contour_show()
fig=plt.figure()
val=cube.stack().max()
plt.imshow(cube.stack(),cmap=cmap,vmin=0,vmax=val)
plt.colorbar()
plt.show()
fig=plt.figure()
plt.imshow(cubelist[2].stack() - cubelist[3].stack(),cmap=cmap,vmin=0,vmax=val)
plt.show()
fig=plt.figure()
for i in range(0,n):
   cb=cubelist[i]
   name=names[i]
   plt.subplot(2,n,i+1)
   plt.title(name+" Fit")
   plt.imshow(cb.stack(),cmap=cmap,vmin=0,vmax=val)
   plt.subplot(2,n,n+i+1)
   plt.title(name+" Residual")
   plt.imshow(res[i],cmap=cmap,vmin=0,vmax=max_re)
plt.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
plt.colorbar(cax=cbar_ax)
plt.show()
plt.pause(1000)
