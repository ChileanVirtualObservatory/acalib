import core.workspace as ws
import numpy as np
import timeit
import cProfile

ws.import_file("fits/M100line.image.fits")
ws.import_file("fits/M100-CO.mom0.fits")
ws.import_file("fits/calibrated.ms.image.spectrum.J113740.6-010454.spw0.image.fits")
ws.import_file("fits/calibrated.ms.line.spw0.source15.image.fits")
ws.import_file("fits/Boom.cm.cln.fits")

elm=ws.elements()
cube=elm['M100line.image-0']
print cube.max()
cube=elm['M100-CO.mom0-0']
print cube.max()
cube=elm['calibrated.ms.image.spectrum.J113740.6-010454.spw0.image-0']
print cube.max()
cube=elm['calibrated.ms.line.spw0.source15.image-0']
print cube.max()
cube=elm['Boom.cm.cln-0']
print cube.max()







