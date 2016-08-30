import sys
sys.path.append('../../')
from acalib import *
from acalib.cupid import pycupid
import astropy.units as u
from astropy.nddata import *
import numpy as np


def struct_builder(caa):
    dims = caa.shape

    #2D data cube
    if len(dims)==2:
        for i in range(dims[0]):
            for j in range(dims[1]):
                pass

    #3D data cube
    elif len(dims)==3:
        for i in range(dims[0]):
            for j in range(dims[1]):
                for z in range(dims[2]):
                    pass


@support_nddata
def clumpfind(data, wcs=None, mask=None, unit=None, rms=0.0):
    data = data.astype(np.float64)
    ret = pycupid.clumpfind(data, dict(), rms)
    return NDData(ret, uncertainty=None, mask=None, wcs=wcs, meta=None, unit=unit)
