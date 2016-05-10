import sys
sys.path.append('../../')

import acalib.io.workspace as ws
import numpy as np
import matplotlib.pyplot as plt
import timeit
import cProfile
import acalib.clumps.wavClumps as wclump
import matplotlib.pyplot as plt
from acalib import acontainer as ac
from acalib.io import graph as gp


cont=ac.AContainer()
cont.load(sys.argv[1])
cube=cont.primary

spar=cube.standarize()

wc = wclump.WavClumps()

#multiresolution data
mr_data = wc.fit(cube)
levels = mr_data.shape[2]

fig = plt.figure()
for i in range(levels):
    plt.imshow(mr_data[:,:,i], origin='image', interpolation="nearest", cmap=plt.cm.gray)
    plt.set_title('Level: {0}'.format(i), fontsize=12)
plt.show()


#newcube=cube.copy()
#newcube.data=caa

#gp.stacked(newcube)
#gp.volume(newcube)
#gp.contour(newcube)

