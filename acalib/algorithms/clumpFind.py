from .. import core
import numpy as np
import pycupid

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
def _clumpfind(data, config, wcs=None, mask=None, unit=None, rms=0.0):
    cube = data
    if len(data.shape) == 4:
        if data.shape[0] == 1:
            cube = data[0,:,:,:]
            if data.shape[1] == 1:
                cube = data[0,:,:]
    elif len(data.shape) == 3:
        if data.shape[0] == 1:
            cube = data[0,:,:]


    ret = pycupid.clumpfind(cube, rms,config=config)
    if ret is not None:
        ret[ret == ret.min()] = 0
        if wcs:
            return NDDataRef(ret, uncertainty=None, mask=None, wcs=wcs, meta=None, unit=unit)
        else:
            return ret
    else:
        return None


class ClumpFind(Algorithm):

    def default_params(self):
        if 'FWHMBEAM' not in self.config:
            self.config['FWHMBEAM'] = 2.0
        if 'VELORES' not in self.config:
            self.config['VELORES'] = 2.0
        if 'ALLOWEDGE' not in self.config:
            self.config['ALLOWEDGE'] = 0
        if 'NAXIS' not in self.config:
            self.config['NAXIS'] = 2
        if 'IDLALG' not in self.config:
            self.config['IDLALG'] = 0
        if 'MINPIX' not in self.config:
            self.config['MINPIX'] = 10


    def run(self, data):
        if type(data) is NDData or type(data) is NDDataRef:
            if len(data.data.shape) > 4:
                raise Exception("Algorithm only support 2D and 3D Matrices")
        else:
            if len(data.shape) > 4:
                raise Exception("Algorithm only support 2D and 3D Matrices")
        # if rms not in config, estimate it
        if 'RMS' not in self.config:
            if type(data) == NDData or type(data)== NDDataRef:
                rms = core.rms(data.data)
            else:
                rms = core.rms(data)
        else:
            rms = self.config['RMS']

        # computing the CAA through clumpfind clumping algorithm
        caa = _clumpfind(data, self.config, rms=rms)

        # computing asocciated structures
        if caa is not None:
            clumps = _struct_builder(caa.data)

            return caa,clumps
        else:
            return None,None
