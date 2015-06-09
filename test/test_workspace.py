import core.workspace as ws
import numpy as np
import timeit
import cProfile

ws.import_file("fits/M100line.image.fits")
ws.import_file("fits/combined-278000.fits")
ws.import_file("fits/calibrated.ms.image.spectrum.J113740.6-010454.spw0.image.fits")
ws.import_file("fits/calibrated.ms.line.spw0.source15.image.fits")
ws.import_file("fits/Boom.cm.cln.fits")

elm=ws.elements()
print " ..."
cube=elm['M100line.image-0']
print " ..."
cube=elm['combined-278000-0']
print " ..."
cube=elm['calibrated.ms.image.spectrum.J113740.6-010454.spw0.image-0']
print " ..."
cube=elm['calibrated.ms.line.spw0.source15.image-0']
print " ..."
cube=elm['Boom.cm.cln-0']







