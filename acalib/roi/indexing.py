import numpy as np
from skimage.filter import threshold_adaptive

from skimage.morphology import binary_opening
from skimage.morphology import disk

from skimage.segmentation import clear_border

from skimage.measure import regionprops
from skimage.measure import label

import matplotlib.pyplot as plt




class SpectraSketcher:
    """
	Create a representation of the cube spectra using pixel samples.

    """

    def __init__(self,nddata):
        """
            Args:
              adata (AData): Datacube to be analysed
        """
        self.cube = nddata

    def cube_spectra(self,samples):
        """
	Create the spectra.
 
        Args:
           samples (int): Number of pixel samples used for the sketch.

        Return:
           spectra (array): An array with the intensity for each frecuency.
           slices  (list):  A list with the slices where emision exist.

        """
        cube = self.cube
        dims = cube.shape
        P_x = dims[2]
        P_x_range = range(P_x)
        P_y = dims[1]
        P_y_range = range(P_y)
        frec = dims[0]

        spectra = np.zeros(frec)

        for i in xrange(samples):
            x_ = np.random.choice(P_x_range,1)
            y_ = np.random.choice(P_y_range,1)
            pixel = cube[:, y_, x_] 
            pixel_masked = self._pixel_processing(pixel)
            spectra += pixel_masked
        spectra = self._pixel_processing(spectra)

        slices = []
        min_slice = -1
        max_slice = -1
        for i in range(frec-1):
            if spectra[i] != 0:
                if min_slice == -1:
                    min_slice = i
                else:
                    if spectra[i+1] == 0:
                        max_slice = i+1
                        slices.append(slice(min_slice,max_slice))
                        min_slice = -1
                    else:
                        if i == frec-2:
                            max_slice = i+1
                            slices.append(slice(min_slice,max_slice))

        return spectra,slices

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


    def vel_stacking(self,data_slice):
        """
            Create an image stacking the frecuency
            
            Args:
           	data_slice: slice object 
            Return:
              image (numpy array): 2D-Array with the stacked cube.

        """
        dims = self.cube.shape
        
	subcube = self.cube[data_slice, :,:]
        stacked = np.sum(subcube,axis=0)
        min_stacked = np.min(stacked)
        
        h = (stacked - min_stacked) / (np.max(stacked) - min_stacked)

        return h


class GaussianSegmentation:
    def __init__(self,prob = 0.05, precision = 2./100):
        self.prob = prob
        self.precision = precision

    def gaussian_mix(self,image):
        prob = self.prob
        dims = image.shape
        rows = dims[0]
        cols = dims[1]
        size = np.min([rows,cols])
        precision = size * self.precision 
        
        image = image.astype('float64')

        w_max = self._optimal_w(image,prob)
        diff = (image - np.min(image)) / (np.max(image)- np.min(image))

        
        g = threshold_adaptive(diff, w_max*w_max,method='mean',offset=0)

        r = w_max/2
        rMin = 2*np.round(precision)


        objects = []
        while (r > rMin):
            background = np.zeros((rows,cols))
            selem = disk(r)
            sub = binary_opening(g,selem)
            sub = clear_border(sub)
            sub = label(sub)
            fts = regionprops(sub)

            plt.imshow(sub, origin='lower')
            plt.show()

            if len(fts) > 0:
                objects.append({'radius': r,'objects':fts})
                for props in fts:
                    C_x, C_y = props.centroid
                    radius = props.equivalent_diameter / 2.
                    kern = 0.01 * np.ones((2*radius, 2*radius))
                    krn = self._kernelsmooth( x = np.ones((2*radius, 2*radius)), kern = kern)
                    krn = np.exp(np.exp(krn))
                    if np.max(krn) > 0:
                        krn = (krn-np.min(krn))/(np.max(krn)-np.min(krn))
                        background = self._kernel_shift(background,krn, C_x, C_y)
            if np.max(background) > 0:
                background = (background-np.min(background))/(np.max(background)-np.min(background))
                diff = diff - background
            diff = (diff-np.min(diff))/(np.max(diff)-np.min(diff))
            g = threshold_adaptive(diff, r*r,method='mean',offset=0)
            r = np.round(r/2.)
        return objects  


    def _optimal_w(self,image,p = 0.05):
        #radiusMin, radius Max and inc in percentages of the image size, p as [0,1] value, image is the original version
        radiusMin = 5
        radiusMax = 40
        inc = 1

        f = (image-np.min(image))/(np.max(image)-np.min(image))
        dims = f.shape
        rows = dims[0]
        cols = dims[1]

        maxsize = np.max([rows,cols])
        imagesize = cols*rows
        radius_thresh = np.round(np.min([rows,cols])/4.)
        unit = np.round(maxsize/100.)

        radiusMin = radiusMin*unit
        radiusMax = radiusMax*unit
        radiusMax = int(np.min([radiusMax,radius_thresh]))
        radius = radiusMin
        inc = inc*unit

        bg = np.percentile(f , p * 100)        
        fg = np.percentile(f, (1-p) * 100)
        min_ov = imagesize

        #print "Radius Min: ", radiusMin, " Radius Max: ", radiusMax, " Increment: ", inc
        #print "Background ", "%.16f "% bg, "Foreground ", "%.16f "%fg
        while(radius <= radiusMax):
            g = threshold_adaptive(f,radius*radius,method='mean',offset=0)
            ov = self._bg_fg(f,g,bg,fg)
            if(ov < min_ov):
                w = radius
                min_ov = ov
            
            radius += inc
        #print "Optimal w ",w
        return w

    def _bg_fg(self,f,g,bg,fg):
        dims = f.shape
        rows = dims[0]
        cols = dims[1]
        fp = 0
        fn = 0
        for rowID in xrange(rows):
            for colId in xrange(cols):
                if g[rowID][colId] == True:
                    if (np.abs(f[rowID][colId]- bg) < np.abs(f[rowID][colId]-fg)):
                        fp += 1
                elif g[rowID][colId] == False:
                    if (np.abs(f[rowID][colId]- bg) > np.abs(f[rowID][colId]-fg)):
                        fn += 1
        overall = fp + fn
        return overall

    def _kernelsmooth(self,x,kern, norm = True):
        # how many rows/cols of zeroes are used to pad.
        width = kern.shape[0]
        pad = int(width/2.)

        #record the width and height the input data matrix
        x_w = x.shape[0]
        x_h = x.shape[1]

        # Are we normalizing the kernel?
        if norm:
            k = kern / np.sum(abs(kern))
        else:
            k = kern

        #Padding with zeros
        x_pad = np.lib.pad(x, ((pad,pad),(pad,pad)), 'constant')

        # Pre-allocate the final (smoothed) data matrix
        s = np.zeros( (x_h, x_w) )

        # Pre-allocate a temporary matrix for the iterative calculations
        temp = np.zeros((width,width))

        # Loop through the data to apply the kernel.
        for col in xrange(x_w):
            for row in xrange(x_h):
                temp = x_pad[row:(row+width),col:(col+width)]
                s[row][col] = np.sum(k*temp)

        return s

    def _kernel_shift(self,back,kernel, x,y):
        rows_back = back.shape[0]
        cols_back = back.shape[1]
        rowsKernel = kernel.shape[0]
        colsKernel = kernel.shape[1]
        rowInit = int(x-rowsKernel/2)
        colInit = int(y-colsKernel/2)

        for row in xrange(rowsKernel):
            for col in xrange(colsKernel):
                if (rowInit + row <= rows_back) and (colInit+col <= cols_back):
                    back[rowInit+row][colInit+col] = kernel[row][col]

        return back
