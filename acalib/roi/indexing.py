import numpy as np


class RoiDetect:
    def __init__(self):
        pass

    def cube_spectra(self,adata,samples):
        cube = adata.data
        dims = adata.shape()
        P_x = dims[2]
        P_x_range = range(P_x)
        P_y = dims[1]
        P_y_range = range(P_y)
        frec = dims[0]


        spectra = np.zeros(frec)

        for i in xrange(samples):
            x_ = np.random.choice(P_x_range,1)
            y_ = np.random.choice(P_y_range,1)
            pixels = cube[:,y_,x_]
            pixel_masked = self._pixel_processing(pixels)
            spectra += pixel_masked
        spectra = self._pixel_processing(spectra)
        return spectra

    def _pixel_processing(self,pixels):
        acum = self._accumulating(pixels)
        diff = self._differenting(acum)
        boxing = self._segmenting(diff)
        boxing = self._erosing(boxing)
        return self._masking(boxing,pixels)

        pass

    def _accumulating(self,pixels):
        return np.cumsum(pixels)
    
    #Can't think in a vectorized way to do it
    def _differenting(self,cumPixels):
        n = len(cumPixels)
        diff = np.zeros(n)
        diff[0] = cumPixels[0]
        for i in xrange(1,n):
            diff[i] = cumPixels[i] - diff[i-1]
        return diff 

    def _segmenting(self,diff):
        n = len(diff)
        boxing = np.zeros(n)
        
        for i in xrange(1,n-1):
            boxing[i] = 1
            if (
                (diff[i] < diff[i-1]) and (diff[i] < diff[i+1])

                or

                (diff[i] > diff[i-1]) and (diff[i] > diff[i+1])
                ):
                boxing[i] = 0

        return boxing

    def _erosing(self, boxing):
        n = len(boxing)
        blocking = np.zeros(n)

        for i in xrange(1,n-1):
            blocking[i] = boxing[i]

            if ( boxing[i-1] == 0 and boxing[i] == 1 and boxing[i+1] == 0 ):
                blocking[i] = 0

        boxing = np.copy(blocking)
        for i in xrange(1, n-1):
            if (blocking[i-1] == 0 and blocking[i]== 1):
                boxing[i-1] = 1
            if (blocking[i]==1 and blocking[i+1]==0):
                boxing[i+1] = 1

        return boxing

    def _masking(self,boxing, pixels):
        n1 = len(boxing)
        n2 = len(pixels)

        if n1 == n2:
            n = n1
            output = np.zeros(n)
            for i in xrange(n):
                output[i] = boxing[i]* pixels[i]
            return output
        else:
            return 0


    def vel_stacking(self,cube,slice_min = None, slice_max = None):
        dims = cube.shape()
        frec = dims[0]
        P_y = dims[1]
        P_x = dims[2]
        if slice_min:
            slice_min = (slice_min, 0, 0)
        else:
            slice_min = (0, 0, 0)

        if slice_max:
            slice_max = (slice_max, P_y, P_x)
        else:
            slice_max = (frec, P_y, P_x)
        
        stacked = cube.stack(lower=slice_min, upper=slice_max)
        min_stacked = min(stacked)
        
        h = (stacked - min_stacked) / (max(stacked) - min_stacked)

        return h