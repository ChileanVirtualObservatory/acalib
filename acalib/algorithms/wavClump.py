import sys
sys.path.append('../../')
import acalib
import astropy.units as u
from astropy.nddata import *
from algorithm import Algorithm
import numpy as np
import numpy.ma as ma
import numba
from joblib import Parallel, delayed


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


def _clumping(data, clumping_method='fellwalker', rms=None):
    if clumping_method=='fellwalker':
        fw = fwalker.FellWalker()
        if rms!=None:
            fw.par['RMS'] = rms
        res = fw.run(data, flat_verification=False)
        del fw
        return res
    elif clumping_method=='clumpfind':
        return None


@numba.jit('boolean (int64[:], int64[:,:])')
def _contained(pixel, clump_pixels):
    retval = False
    for pp in clump_pixels:
        if (pp[0]==pixel[0] and pp[1]==pixel[1] and pp[2]==pixel[2]):
            retval = True
            break
    return retval


def _empty_tree(indexes):
    tree = dict()
    for ind in indexes:
        tree[ind] = list()
    return tree




class WavClump(Algorithm):
    
    def default_params(self):
        if 'wavelet' not in self.config:
            self.config['wavelet'] = 'sym5'
        if 'nlevel' not in self.config:
            self.config['nlevel'] = 4
        if 'clumping_method' not in self.config:
            self.config['clumping_method'] = 'fellwalker'
        if 'keep_rms' not in self.config:
            self.config['keep_rms'] = False
        if 'n_jobs' not in self.config:
            self.config['n_jobs'] = 2

    
    def mra(self, ml_data):
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
    
    
    def run(self, np_data, ml_data):
        #perform multiresolution analysis
        np_approx, ml_approx = self.mra(ml_data)
        np_approx.insert(0, np_data)
        ml_approx.insert(0, ml_data)
        
        #compute RMS of original data
        if self.keep_rms:
            rms = np.sqrt(np.sum(np_data*np_data) / np.size(np_data))
        
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