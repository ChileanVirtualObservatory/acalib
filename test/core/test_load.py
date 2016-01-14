import core.workspace as ws
import numpy as np
import matplotlib.pyplot as plt
import timeit
import cProfile
import clumps.fellWalker as fwalker
import matplotlib.pyplot as plt

ws.import_file("fits/M100line.image.fits")
elm=ws.elements()
cube=elm['M100line.image-0']
print "done!"
