import numpy as np
from acalib.cupid import pycupid
from astropy.nddata import *

from .algorithm import Algorithm


# storing unusable pixels for now (-1)
def _struct_builder(caa):
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
def _fellwalker(data, config, wcs=None, mask=None, unit=None, rms=0.0):
    data = data.astype(np.float64)
    ret = pycupid.fellwalker(data, config, rms)
    return NDData(ret, uncertainty=None, mask=None, wcs=wcs, meta=None, unit=unit)


class FellWalker(Algorithm):

    def default_params(self):
        if 'FWHMBEAM' not in self.config:
            self.config['FWHMBEAM'] = 2.0
        if 'VELORES' not in self.config:
            self.config['VELORES'] =  2.0  

    def run(self, data):
        # if rms not in config, estimate it
        if 'RMS' not in self.config:
            rms = rms(data)
        else:
            rms = self.config['RMS']

        # computing the CAA through CUPID's fellwalker clumping algorithm
        caa = _fellwalker(data, self.config, rms=rms)

        # computing asocciated structures
        clumps = _struct_builder(caa.data)

        return caa,clumps
