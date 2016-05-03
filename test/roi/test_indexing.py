import sys
sys.path.append('../../')

import acalib
import acalib.io.formats as io
import acalib.roi as ri

c = acalib.AContainer()
io.load_to_cont('/home/cvalenzu/Escritorio/indexing_utfsm/ALMA00000085.fits',c)
cube = c.primary
rd = ri.RoiDetect()
spectra = rd.cube_spectra(cube,200)

h1 = rd.vel_stacking(cube,67,96)