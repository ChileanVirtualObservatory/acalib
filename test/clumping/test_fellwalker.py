import sys
sys.path.append('../../')


import acalib.vo.workspace as ws
import numpy as np
import matplotlib.pyplot as plt
import timeit
import cProfile
import acalib.clumps.fellWalker as fwalker
import matplotlib.pyplot as plt

binpath='../../bindata/fits/cubes/'
ws.import_file(binpath+"M100line.image.fits")
#ws.import_file(binpath+"Orion.methanol.cbc.contsub.image.fits")
#ws.import_file(binpath+"Boom.cm.cln.fits")
#ws.import_file(binpath+"Antennae_North.CO3_2Line.Clean.pcal1.image.fits")
#ws.import_file(binpath+"Antennae_South.CO3_2Line.Clean.pcal1.image.fits")
#ws.import_file(binpath+"calibrated.ms.contsub.bin4.line.fits")
#ws.import_file(binpath+"calibrated.ms.contsub.bin4.line.image.fits")

elm=ws.elements()
cube=elm['M100line.image-0']
#cube=elm['Orion.methanol.cbc.contsub.image-0']
#cube=elm['Boom.cm.cln-0']
#cube=elm['Antennae_North.CO3_2Line.Clean.pcal1.image-0']
#cube=elm['Antennae_South.CO3_2Line.Clean.pcal1.image-0']
#cube=elm['calibrated.ms.contsub.bin4.line-0']
#cube=elm['calibrated.ms.contsub.bin4.line.image-0']

spar=cube.standarize()

fw = fwalker.FellWalker()

caa,clump=fw.fit(cube)

newcube=cube.copy()
newcube.data=caa
newcube.stacked_show()
newcube.volume_show()
newcube.contour_show()

