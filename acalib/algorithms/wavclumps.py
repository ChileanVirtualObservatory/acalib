import numba
import numpy as np
import numpy.ma as ma
from astropy.nddata import *
from joblib import Parallel, delayed
from pycupid import clumpfind, fellwalker
from .. import core


def _struct_builder(data, caa):
    """
    TODO: docstring
    """
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

    peaks = {}
    for i,pl in clumps.items():
        peaks[i] = None
        max_value = -float("inf")
        # for pixel position in pixel list
        for pp in pl:
            if data[pp] > max_value:
                max_value = data[pp]
                peaks[i] = pp

    return clumps,peaks


def _clumping(data, clumping_method='fellwalker', clumping_conf=dict(), rms=None):
    """
    TODO: docstring
    """
    if rms is None:
        rms = acalib.rms(data)

    if clumping_method=='fellwalker':
        caa = fellwalker(data, rms, config=clumping_conf)
    elif clumping_method=='clumpfind':
        caa = clumpfind(data, rms, config=clumping_conf)
    else:
        return None

    clumps,peaks = _struct_builder(data, caa)
    return caa,clumps,peaks,rms


@numba.jit('boolean (int64[:], int64[:,:])')
def _contained(pixel, clump_pixels):
    """
    TODO: docstring
    """
    retval = False
    for pp in clump_pixels:
        if (pp[0]==pixel[0] and pp[1]==pixel[1] and pp[2]==pixel[2]):
            retval = True
            break
    return retval


def _empty_tree(indexes):
    """
    TODO: docstring
    """
    tree = dict()
    for ind in indexes:
        tree[ind] = list()
    return tree



class WavClumps(Algorithm):
    
    def default_params(self):
        if 'WAVELET' not in self.config:
            self.config['WAVELET'] = 'sym5'
        if 'NLEVEL' not in self.config:
            self.config['NLEVEL'] = 4
        if 'CLUMPING_METHOD' not in self.config:
            self.config['CLUMPING_METHOD'] = 'fellwalker'
        if 'KEEP_RMS' not in self.config:
            self.config['KEEP_RMS'] = False
        if 'NJOBS' not in self.config:
            self.config['NJOBS'] = 2

    
    def mra(self, ml_data, eng):
        """
        Computes approximations and details for each level
        """
        ml_approx = eng.mra(ml_data, self.wavelet, self.nlevel, nargout=1)
        
        #casting to numpy.ndarray
        np_approx = map(np.asarray, ml_approx)
        for i in range(len(np_approx)):
            np_approx[i] = ma.masked_array(np_approx[i], mask=np.isnan(np_approx[i]))
        return (np_approx, ml_approx)
        
    
    def tree_build(self, caa, clumps, peaks):
        #build the hierarchical tree
        htree = []
        for lvl in range(self.nlevel-1,-1,-1):
            tree = _empty_tree(clumps[lvl+1].keys())
            for ind1,pk in peaks[lvl].items():
                for ind2,clump in clumps[lvl+1].items():
                    if _contained(pk,clump):
                        tree[ind2].append(ind1)
            htree.append(tree)
        return htree
    
    
    def run(self, np_data, ml_data, eng):
        #perform multiresolution analysis
        np_approx, ml_approx = self.mra(ml_data, eng)
        np_approx.insert(0, np_data)
        ml_approx.insert(0, ml_data)
        
        #compute RMS of original data
        if self.keep_rms:
            rms = acalib.rms(np_data)
        
        #perform clumping algorithm over original data and each approximation level
        if self.keep_rms:
            parallelizer = Parallel(n_jobs=self.n_jobs)
            # this iterator returns the functions to execute for each task
            tasks_iterator = [delayed(_clumping)(data, rms=rms) for data in np_approx]
            result = parallelizer(tasks_iterator)            
        else:
            parallelizer = Parallel(n_jobs=self.n_jobs)
            # this iterator returns the functions to execute for each task
            tasks_iterator = [delayed(_clumping)(data) for data in np_approx]
            result = parallelizer(tasks_iterator)
        
        #ordering results
        caa = list()
        clumps = list()
        peaks = list()
        rmss = list()
        
        for res in result:
            caa.append(res[0])
            #parsing clumps pixel lists to numpy arrays
            for index in res[1]:
                pixel_list = res[1].pop(index)
                res[1][index] = np.asarray(pixel_list)
            clumps.append(res[1])
            #parsing peaks positions (tuples) to numpy arrays
            for index in res[2]:
                ppeak = res[2].pop(index)
                res[2][index] = np.asarray(ppeak)
            peaks.append(res[2])
            rmss.append(res[3])
        
        #build the hierarchical tree
        htree = self.tree_build(caa,clumps,peaks)
        
        max_pixels = []
        mean_pixels = []
        #computing clump associated structures
        for i in range(self.nlevel+1):
            max_val = 0
            acc_val = 0
            for pixels in clumps[i].values():
                npix = len(pixels)
                acc_val += npix
                if npix>max_val: max_val = npix
            acc_val /= len(clumps[i])
            max_pixels.append(max_val)
            mean_pixels.append(acc_val)
            
        #some aditional information of each approximation
        entropy = []
        variance = []
        for i in range(self.nlevel+1):
            entropy.append(eng.entropy(ml_approx[i], nargout=1))
            variance.append(np.std(np_approx[i])**2)
        
        #storing results
        self.caa = caa
        self.clumps = clumps
        self.peaks = peaks
        self.rmss = rmss
        self.htree = htree
        self.max_pixels = max_pixels
        self.mean_pixels = mean_pixels
        self.entropy = entropy
        self.variance = variance
        self.mra_data = np_approx
        return caa,clumps
    
    def summarize(self):
        for i in range(self.nlevel+1):
            print('Level {0}:'.format(i))
            print('Number of clumps: {0}'.format(len(self.clumps[i])))
            print('RMS of data cube: {0}'.format(self.rmss[i]))
            print('Number of pixels of bigger clump: {0}'.format(self.max_pixels[i]))
            print('Mean number of pixels of clums: {0}'.format(self.mean_pixels[i]))
            print('Entropy of approximation: {0}'.format(self.entropy[i]))
            print('Variance of approximation: {0}'.format(self.variance[i]))
            print('------------------------------------------')