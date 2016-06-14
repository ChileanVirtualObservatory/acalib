import sys
sys.path.append('../../')

import acalib.io.workspace as ws
import numpy as np
import matplotlib.pyplot as plt
import timeit
import cProfile
import acalib.clumps.clumpFind as cfind
import matplotlib.pyplot as plt
from acalib import acontainer as ac
#from acalib.io import graph as gp


cont=ac.AContainer()
cont.load(sys.argv[1])
cube=cont.primary

spar=cube.standarize()

cf = cfind.ClumpFind()

caa,clumps = cf.fit(cube)

#newcube=cube.copy()
#newcube.data=caa

#gp.stacked(newcube)
#gp.volume(newcube)
#gp.contour(newcube)

