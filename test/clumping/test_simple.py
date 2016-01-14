import core.workspace as ws
import numpy as np
import matplotlib.pyplot as plt
import timeit
import cProfile
import clumps.simpleClumps as sclumps
import matplotlib.pyplot as plt

ws.import_file("fits/M100line.image.fits")

elm=ws.elements()
cube=elm['M100line.image-0']
spar=cube.standarize()
sc=sclumps.SimpleClumps()
print "spar",spar

# Fitme :S
clist=sc.fit(cube,verbose=True)
plt.subplot(1, 3, 1)
plt.imshow(cube.get_stacked())
plt.subplot(1, 3, 2)
plt.imshow(sc.data.get_stacked())
plt.subplot(1, 3, 3)
plt.imshow(sc.syn.get_stacked())
plt.show()




