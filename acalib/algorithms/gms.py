import acalib

import numpy as np
try:
    from skimage.filters import threshold_local
except:
    from skimage.filter import threshold_local
from skimage.measure import label,regionprops
from skimage.morphology import binary_opening, disk
from skimage.segmentation import clear_border


from astropy.nddata import support_nddata
from .algorithm import Algorithm

from acalib.core.analysis import _optimal_w, _kernelsmooth, _kernel_shift


@support_nddata
def get_data(data,wcs = None):
    return data,wcs

class GMS(Algorithm):
    """
    Gaussian Multiscale Segmentation:

    Using a mixture of gaussians make an multiscale segmentation to get the region of interest of a 2D astronomical image.
    
    Parameters
    ----------
    params : dict (default = None)
        Algorithm parameters, allowed keys:
         
        P : float (default = 0.05)
            Thresholding quantile for multiscale segmentation.
        PRECISION : float (default = 0.02)
            Smallest scale percentage for the multiscale segmentation.

    References
    ----------
    .. [1] Araya, M., Candia, G., Gregorio, R., Mendoza, M., & Solar, M. (2016). Indexing data cubes for content-based searches in radio astronomy. Astronomy and Computing, 14, 23-34.   
    """
    def default_params(self):
        if 'P' not in self.config:
            self.config['P'] = 0.05
        if 'PRECISION' not in self.config:
            self.config['PRECISION'] = 0.02


    def run(self, data):
        """        
        Parameters
        ----------
        data : (M,N) numpy.ndarray or astropy.nddata.NDData or astropy.nddata.NDDataRef
            Velocity collapsed image
        
        Returns
        ----------
        List of labeled images. 
        """
        data,wcs = get_data(data)

        #TODO: check for wcs != None
        if len(data.shape) > 2:
            log.error("Only 2D images supported")
            raise ValueError("Only 2D images supported")

        image_list = []

        image = data
        image[np.isnan(image)] = 0

        prob = self.config["P"]
        dims = image.shape
        rows = dims[0]
        cols = dims[1]
        size = np.min([rows, cols])
        precision = size * self.config["PRECISION"]

        image = image.astype('float64')

        #Getting optimal radius for first step segmentation
        w_max = _optimal_w(image, prob)

        diff = (image - np.min(image)) / (np.max(image) - np.min(image))

        tt = w_max * w_max

        # Initial segmentation
        if tt % 2 == 0:
            tt += 1
        adaptive_threshold = threshold_local(diff, int(tt), method='mean', offset=0)#(diff, int(tt), offset=0)
        g = diff > adaptive_threshold

        r = w_max / 2

        #Smallest radious for region
        rMin = 2 * np.round(precision)

        #Iterate over radious dividing it by 2 
        # in each iteration
        while (r > rMin):
            background = np.zeros((rows, cols))
            
            #Label the previous segmentation
            #calculate shape features for each 
            #connected region.
            selem = disk(r)
            sub = binary_opening(g, selem)
            sub = clear_border(sub)
            sub = label(sub)
            fts = regionprops(sub)

            #image_list.append(NDData(sub, wcs=wcs))
            # Non NNData version (without wcs... lets check if it pass)
            image_list.append(sub)

            #Uses a gaussian mix to fit the region 
            #then smooth it and remove the gaussian mixture
            #from the image, uses this new image to continue in next 
            #iteration
            if len(fts) > 0:
                for props in fts:
                    C_x, C_y = props.centroid

                    radius = int(props.equivalent_diameter / 2.)
                    kern = 0.01 * np.ones((2 * radius, 2 * radius))
                    krn = _kernelsmooth(x=np.ones((2 * radius, 2 * radius)), kern=kern)
                    krn = np.exp(np.exp(krn))
                    if np.max(krn) > 0:
                        krn = (krn - np.min(krn)) / (np.max(krn) - np.min(krn))
                        background = _kernel_shift(background, krn, C_x, C_y)
            if np.max(background) > 0:
                background = (background - np.min(background)) / (np.max(background) - np.min(background))
                diff = diff - background
            diff = (diff - np.min(diff)) / (np.max(diff) - np.min(diff))
            tt = int(r * r)
            if tt % 2 == 0:
                tt += 1
            adaptive_threshold = threshold_local(diff, tt, offset=0, method='mean')
            g = diff > adaptive_threshold

            r = np.round(r / 2.)

        return image_list


