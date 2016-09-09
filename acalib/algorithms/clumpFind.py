1import sys
sys.path.append('../../')
from acalib import *
from acalib.cupid import pycupid
import astropy.units as u
from astropy.nddata import *
import numpy as np


def struct_builder(caa):
    dims = caa.shape
    clumps = dict()

    #2D data cube
    if len(dims)==2:
        for i in range(dims[0]):
            for j in range(dims[1]):
                if caa[i,j] in clumps:
                    clumps[caa[i,j]].append((i,j))
                else:
                    clumps[caa[i,j]] = [(i,j)]
    #3D data cube
    elif len(dims)==3:
        for i in range(dims[0]):
            for j in range(dims[1]):
                for k in range(dims[2]):
                if caa[i,j,k] in clumps:
                    clumps[caa[i,j,k]].append((i,j,k))
                else:
                    clumps[caa[i,j,k]] = [(i,j,k)]
    return clumps


@support_nddata
def _clumpfind(data, wcs=None, mask=None, unit=None, rms=0.0):
    data = data.astype(np.float64)
    ret = pycupid.clumpfind(data, dict(), rms)
    return NDData(ret, uncertainty=None, mask=None, wcs=wcs, meta=None, unit=unit)


class ClumpFind:

    def __init__(self, config):
        if config is not None:
            self.default_params()

    def default_params(self, config):
        pass

    def run(self, data):
        ret = _clumpfind(data, rms)
