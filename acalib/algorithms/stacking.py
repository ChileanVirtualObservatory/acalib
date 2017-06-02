import acalib
from .algorithm import Algorithm

from numpy import mean
from scipy.stats import signaltonoise
from astropy.nddata import NDData,NDDataRef

class Stacking(Algorithm):
    """
    Create a stacked image using a template image and a set of different images from same object.
    """

    def default_params(self):
        pass

    def run(self, template_data, data_cont):
        """
            Run the stacking algorithm given a template image and a container of images.

            Parameters
            ----------
            template_data : (M,N) numpy.ndarray
                Astronomical image.
            data_cont : acalib.container
                An images set container

            Returns
            -------
            result : (M,N) numpy.ndarray
                Image stacked
        """
        if template_data is NDData or template_data is NDDataRef:
            template_data = template_data.data

        for i in range(len(data_cont.images)):
            if data_cont.images[i] is not NDData or data_cont.images[i] is not NDDataRef:
                data_cont.images[i] = NDDataRef(data_cont.images[i])

        tprops = acalib.core.transform.fits_props(template_data)
        scaled = acalib.core.transform.scale(data_cont, tprops['major'])
        rotated, angles = acalib.core.transform.rotate(scaled, tprops['angle'])
        aligned = acalib.core.transform.crop_and_align(rotated, angles)
        result = mean(aligned , axis=0)

        
        return result