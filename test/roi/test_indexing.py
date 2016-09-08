import sys
sys.path.append('../../')

import acalib
from astropy.nddata import NDData


c = acalib.Container()
c.load_fits('../../../bindata/fits/cubes/ALMA01000740.fits')
cube = c.primary

spectra,slices = acalib.cube_spectra(cube,5000)

sl = slice(26, 48)
image = acalib.vel_stacking(cube,sl)

objects,images = acalib.gaussian_mix(image, images=True)