import sys
sys.path.append('../../')
from acalib import *
from acalib.cupid import pycupid
import astropy.units as u
from astropy.nddata import *

@support_nddata
def clumpfind(data, wcs=None, mask=None, unit=None, rms=0.0):
    ret = pycupid.clumpfind(data, dict(), rms)
    return NDData(ret, uncertainty=None, mask=None, wcs=wcs, meta=None, unit=unit)