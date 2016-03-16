import sys
sys.path.append('../../')


import acalib.io.workspace as ws
import numpy as np
import timeit
import cProfile

binpath='../../bindata/fits/cubes/'
ws.import_file(binpath+"M100line.image.fits")
ws.import_file(binpath+"Orion.methanol.cbc.contsub.image.fits")
ws.import_file(binpath+"Boom.cm.cln.fits")

elm=ws.elements()
cube=elm['M100line.image-0']
print cube.max()
cube=elm['Orion.methanol.cbc.contsub.image-0']
print cube.max()
cube=elm['Boom.cm.cln-0']
print cube.max()







