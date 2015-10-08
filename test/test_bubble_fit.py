import core.workspace as ws
import time
import numpy as np
import matplotlib.pyplot as plt
import timeit
import cProfile
import clumps.bubbleClumps as cl
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

porc=[]
elem=[]
snr=[]
tim=[]


total=cube.sum()
telem=float(cube.count())
n=6
for i in range(n):
   snrlimit=2.0 - 0.35*i
   snr.append(snrlimit)
   gc.par['SNRLIMIT']=snrlimit
   start = time.clock()
   gc.fit(cube,verbose=True)
   end = time.clock()
   elem.append(len(gc.amplitudes)/telem)
   porc.append(gc.syn.sum()/total)
   tim.append(end-start)
   plt.subplot(2, n,0*n + i + 1)
   plt.imshow(gc.syn.get_stacked())
   plt.subplot(2, n,1*n + i + 1)
   plt.imshow(gc.data.get_stacked())
   #plt.subplot(3, n,2*n + i + 1)
   #pos=np.array(gc.positions)
   #plt.scatter(pos[:,2], pos[:,1],marker='x')

print "STATS"
print "elem",elem
print "porc",porc
print "time",tim

#gc.par['NRMS']=1.0
#gc.fit(cube,verbose=True)

plt.show()

plt.clf()
plt.subplot(1, 3,1)
plt.plot(snr,elem,'k')
plt.subplot(1, 3,2)
plt.plot(snr,porc,'b')
plt.subplot(1, 3,3)
plt.plot(snr,tim,'r')


plt.show()

