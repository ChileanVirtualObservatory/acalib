import acalib
from .algorithm import Algorithm

from numpy import mean
from scipy.stats import signaltonoise


class Stacking(Algorithm):
    def default_params(self):
    	pass

    def run(self, template_data, data_cont):

        tprops = acalib.core.transform.fits_props(template_data)
        scaled = acalib.core.transform.scale(data_cont, tprops['major'])
        rotated, angles = acalib.core.transform.rotate(scaled, tprops['angle'])
        aligned = acalib.core.transform.crop_and_align(rotated, angles)
        result = mean(aligned , axis=0)
        return result, signaltonoise(result)