import sys
import pywt
import numpy as np
import matplotlib.pyplot as plt


class WavClumps:
    def __init__(self):
        self.defaultParams()
    
    def defaultParams(self):
        self.par = dict()
        
        """ Generic parameters  """
        # Spectral resoluion in pixels
        self.par['VELORES'] = 2.0
        # Beam resolution in pixels
        self.par['FWHMBEAM'] = 2.0
        # Maximum Clumps
        self.par['MAXCLUMPS'] = sys.maxint
        # The lower threshold for clump values to a user-specified multiple of the RMS noise.
        self.par['THRESH'] = 1.0

    #performs multiresolution analysis through Stationary Wavelet Tranform
    def mr_analysis(self, data, wavelet='bior1.3', max_level=None):
        (M,N) = data.shape
        _max_level = pywt.swt_max_level(min(M,N))
        #setting max allowed level
        if max_level==None or max_level>_max_level:
            max_level = _max_level

        #performing the SWT at different scales
        mr_swt = pywt.swt2(data, wavelet, level=max_level, start_level=0)

        #allocating memory for results (image at different resolutions) and storing it
        res = np.empty((M,N,max_level))
        level_ind = 0
        for LL,_ in mr_swt:
            res[:,:,level_ind] = LL
            level_ind += 1
        return res


    def fit(self, cube):
        #storing original cube
        orig_cube = cube.copy()

        # Set the RMS, or automatically find an estimate for it
        if not self.par.has_key('RMS'):
            rms = cube.estimate_rms()
            self.par['RMS'] = rms

        #stacking cube on Velocity axis
        data = cube.stack(axis=(0))

        #performing multiresolution analysis
        mr_data = mr_analysis(data, max_level=10)

        return mr_data

        


