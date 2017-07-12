import acalib
from .algorithm import Algorithm

from numpy import mean
from scipy.stats import signaltonoise
import scipy.ndimage as scnd
from astropy.nddata import NDData,NDDataRef

class Stacking(Algorithm):
    """
    Create a stacked image using a template image and a set of different images from same object.
    """

    def default_params(self):
        pass

    def run(self, template_data, images):
        """
            Run the stacking algorithm given a template image and a container of images.

            Parameters
            ----------
            template_data : (M,N) numpy.ndarray
                Astronomical image.
            images : list of (M,N) numpy.ndarray
                A list of images.

            Returns
            -------
            result : (M,N) numpy.ndarray
                Image stacked
        """
        if type(template_data) is NDData or type(template_data) is NDDataRef:
            template_data = template_data.data

        for i in range(len(images)):
            if type(images[i]) is not NDData or type(images[i]) is not NDDataRef:
                images[i] = NDDataRef(images[i])

        tprops = acalib.core.transform.fits_props(template_data)

        # TODO: Replace with core.transform.scale once it stops using
        # a acalib.container.
        majorAxisTemplate = tprops['major']
        scaledData = []
        for i in range(len(images)):
            prop = acalib.core.transform.fits_props(images[i].data)
            sc = majorAxisTemplate / prop['major']
            scaledData.append(scnd.zoom(prop['orig'], sc))
        scaled = scaledData

        rotated, angles = acalib.core.transform.rotate(scaled, tprops['angle'])
        aligned = acalib.core.transform.crop_and_align(rotated, angles)
        result = mean(aligned,axis=0)

        return result
