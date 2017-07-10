from acalib import upi
from .algorithm import Algorithm

from acalib import core
import numpy as np
from astropy.table import Table
from astropy import log

def _vertical_flux_decomposition(rep,delta,noise,kernel,n_partitions,shape):
    n_rep=len(rep)/n_partitions
    img_list=[]
    vmax=0
    for i in range(n_partitions):
        synNew=np.zeros(shape)
        ini=n_rep*i
        end=n_rep*(i+1) - 1
        p_rep=rep[ini:end]
        core.synthesize_bubbles(synNew,p_rep,kernel,noise,delta)
        img=synNew.sum(axis=(0))
        vmax=max(vmax,img.max())
        img_list.append(img)
    return img_list,vmax



class HRep(Algorithm):
    """
    Create a stacked image using a template image and a set of different images from same object.
    """

    def toVFD(table, n_partitions, shape):
        rep, delta, noise, kernel = HRep.toTuple(table)
        return _vertical_flux_decomposition(rep, delta, noise, kernel, n_partitions, shape)

    def toImage(table, template):
        rep, delta, noise, kernel = HRep.toTuple(table)
        synNew = np.zeros(template.data.shape)
        core.synthesize_bubbles(synNew, rep, kernel, noise, delta)
        scale = table.meta['SCALE']
        shift = table.meta['SHIFT']
        return upi.Data(synNew * scale - shift, wcs=template.wcs, unit=template.unit, meta=template.meta)

    def toTuple(table):
        rep = np.array([table[c] for c in table.columns])
        rep = rep.T
        delta = np.array([table.meta['DELTAX'], table.meta['DELTAY'], table.meta['DELTAZ']])
        noise = table.meta['NOISE']
        P = core.precision_from_delta(delta, 0.1)
        kernel = core.create_mould(P, delta)
        return rep, delta, noise, kernel

    def default_params(self):
        if 'KERNEL' not in self.config:
            self.config['KERNEL'] = 'PIXEL'
        if 'DELTA' not in self.config:
            self.config['DELTA'] = None
        if 'NOISE' not in self.config:
            self.config['RMS'] = None
        if 'SNR' not in self.config:
            self.config['SNR'] = None
        if 'STANDARIZE' not in self.config:
            self.config['STANDARIZE'] = True
        if 'VERBOSE' not in self.config:
            self.config['VERBOSE'] = False
        if 'GAMMA' not in self.config:
            self.config['GAMMA'] = 0.1

    def run(self, cube):
        """
            Run the Homogenous Representation algorithm a Data Object.

            Parameters
            ----------
            cube : the cube to represent

            Returns
            -------
            result : ???
        """
        delta = self.config['DELTA']
        noise = self.config['NOISE']
        snr = self.config['SNR']
        standarize = self.config['STANDARIZE']
        verbose = self.config['VERBOSE']
        gamma = self.config['GAMMA']

        scale = 1.0
        shift = 0.0

        if standarize:
            (cube, scale, shift) = upi.standarize(cube)

        if snr is None:
            snr = core.snr_estimation(cube.data, mask=cube.mask, noise=noise)

        if delta is None:
            if cube.meta is None:
                delta = [1, 1, 1]
            else:
                spa = np.ceil((np.abs(cube.meta['BMIN'] / cube.meta['CDELT1']) - 1) / 2.0)
                delta = [1, spa, spa]
        if noise is None:
            noise = core.rms(cube.data)

        if self.config['KERNEL'] == 'PIXEL':
            positions,synthetic,residual=core.scat_pix_detect(cube.data,threshold=snr*noise,full_output=True)

        if self.config['KERNEL'] == 'METABUBBLE':

            # if verbose:
            #     log.info(snr, noise, delta)
            P = core.precision_from_delta(delta, gamma)
            kernel = core.create_mould(P, delta)
            sym = core.eighth_mould(P, delta)
            positions, synthetic, residual, energy, elist = core.scat_kernel_detect(cube.data,delta=delta,kernel=kernel,threshold=snr*noise,noise=noise,full_output=True,sym=sym,verbose=verbose)
        positions = np.array(positions)
        # Pack metadata
        metapack = dict()
        metapack['DELTAX']=delta[0]
        metapack['DELTAY'] = delta[1]
        metapack['DELTAZ'] = delta[2]
        metapack['SCALE'] = scale
        metapack['SHIFT'] = shift
        metapack['KERNEL'] = self.config['KERNEL']
        metapack['NOISE'] = noise
        metapack['SNR'] = snr
        metapack['GAMMA'] = gamma
        rep = Table(positions, names=['x','y','z'],meta=metapack)

        return rep, upi.Data(synthetic,meta=cube.meta,mask=cube.mask,unit=cube.unit,wcs=cube.wcs),upi.Data(residual,meta=cube.meta,mask=cube.mask,unit=cube.unit,wcs=cube.wcs)
