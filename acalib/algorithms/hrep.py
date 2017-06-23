import acalib
from .algorithm import Algorithm
from astropy.nddata import support_nddata

# TODO: This is non-generic. Try to use the UPI (it can be done!)
@support_nddata
def vel_stacking(data,data_slice,wcs=None,uncertainty=None, mask=None, meta=None, unit=None):
     """
     Create an image collapsing the frecuency axis

     Parameters
     ----------
     data : numpy.ndarray or astropy.nddata.NDData or astropy.nddata.NDDataRef
         Astronomical 2D image

     slice : slice object
         Sector to be collapsed

     Returns
     -------
     image (NDDataRef): 2D-Array with the stacked cube.

     """
     if len(data.shape) != 3:
         log.error("Cube needs to be a 3D array")
         raise ValueError("Cube needs to be a 3D array")
     dims = data.shape
     subcube = data[data_slice, :,:]
     stacked = np.sum(subcube,axis=0)
     if wcs:
         wcs = wcs.dropaxis(2)

         return NDDataRef(stacked, uncertainty=uncertainty, mask=mask,wcs=wcs, meta=meta, unit=unit)
     else:
         return stacked


class ScatterPixelRepresentation(Algorithm):

    def default_params(self):
        if 'P' not in self.config:
            self.config['P'] = 0.05
        if 'PRECISION' not in self.config:
            self.config['PRECISION'] = 0.02
        if 'RANDOM_STATE' not in self.config:
            self.config['RANDOM_STATE'] = None
        if 'SAMPLES' not in self.config:
            self.config["SAMPLES"] = 1000

    def run(self, data):

        if data.wcs:
            wcs = data.wcs
        else:
            wcs = None

        c = acalib.Container()

        spectra, slices = acalib.core.spectra_sketch(data.data, self.config["SAMPLES"], self.config["RANDOM_STATE"])

        pp_slices = []
        for slice in slices:
            pp_slice = vel_stacking(data.data, slice)
            labeled_images = acalib.core.gaussian_mix(pp_slice, prob=self.config["P"],
                                                      precision=self.config["PRECISION"])

            if wcs is not None:
                freq_min = float(wcs.all_pix2world(0, 0, slice.start, 1)[2])
                freq_max = float(wcs.all_pix2world(0, 0, slice.stop, 1)[2])
            else:
                freq_min = None
                freq_max = None

            table = acalib.core.measure_shape(pp_slice, labeled_images, freq_min, freq_max)
            if len(table) > 0:
                c.tables.append(table)
                print("pp")
                c.images.append(pp_slice)
                print(type(pp_slice))
                print("lb")
                c.images.extend(labeled_images)
                print(type(labeled_images))
        c.images.insert(0, data)
        c.primary = c.images[0]
        return c
