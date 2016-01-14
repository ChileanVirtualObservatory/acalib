import core.workspace as ws
import numpy as np
import matplotlib.pyplot as plt
import timeit
import cProfile

ws.import_file("fits/M100line.image.fits")

elm=ws.elements()
cube=elm['M100line.image-0']
y,index=cube.max()
print "max:",y, index
pos=cube.index_to_wcs(index)
print "wcspos:", pos 
low,up=cube.index_from_window(pos,[0.01,0.01,20000000])
print "low/up:",low,up
print "features"
adn=cube.get_features(low,up)
print adn
yyy=cube.get_slice(low,up).ravel()
print "vals"
print yyy
plt.imshow(cube.get_stacked())
plt.show()





