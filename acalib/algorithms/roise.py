import acalib
from .algorithm import Algorithm
from .gms import GMS
from astropy.nddata import support_nddata, NDDataRef, NDData
import numpy as np
from acalib.upi import Data
from astropy.table import Table

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

         return Data(stacked, uncertainty=uncertainty, mask=mask,wcs=wcs, meta=meta, unit=unit)
     else:
         return stacked


@support_nddata
def measure_shape(data, labeled_images, min_freq=None, max_freq=None, wcs=None):
    """ Measure a few statistics from labeled images """
    # TODO: Document this function
    objects = list()
    intensity_image = data
    for image in labeled_images:
        objs_properties = acalib.core.get_shape(image, intensity_image)
        objects.extend(objs_properties)

    if len(objects) == 0:
        return Table()

    names = ["CentroidRa", "CentroidDec", "MajorAxisLength", "MinorAxisLength",
             "Area", "Eccentricity", "Solidity", "FilledPercentaje", "MaxIntensity", "MinIntensity", "AverageIntensity"]

    meta = {"name": "Object Shapes"}

    if min_freq is not None:
        meta["min_freq_hz"] = min_freq

    if max_freq is not None:
        meta["max_freq_hz"] = max_freq

    t = Table(rows=objects, names=names, meta=meta)
    return t

class RoiSE(Algorithm):
    """
    Perform an unsupervised region of interest detection and extract shape features.

    Parameters
    ----------
    params : dict (default = None)
        Algorithm parameters, allowed keys:    
        
        P : float (default = 0.05)
            Thresholding quantile for multiscale segmentation.
        PRECISION : float (default = 0.02)
            Smallest scale percentage for the multiscale segmentation.
        SAMPLES : int (default = 1000)
            Number of pixels used to generate the spectra sketch.
        RANDOM_STATE : int (default = None)
            Seed for random smpling. 


    References
    ----------
    
    .. [1] Araya, M., Candia, G., Gregorio, R., Mendoza, M., & Solar, M. (2016). Indexing data cubes for content-based searches in radio astronomy. Astronomy and Computing, 14, 23-34.
    
    """
    def default_params(self):
        if 'P' not in self.config:
            self.config['P'] = 0.05
        if 'PRECISION' not in self.config:
            self.config['PRECISION'] = 0.02
        if 'RANDOM_STATE' not in self.config:
            self.config['RANDOM_STATE'] = None
        if 'SAMPLES' not in self.config:
            self.config["SAMPLES"] = 1000


    def run(self, cube):
        """
            Run the indexing algorithm on a given data cube.

            Parameters
            ----------
            data : (M,N,Z) numpy.ndarray or astropy.nddata.NDData or astropy.nddata.NDDataRef
                Astronomical data cube.

            Returns
            -------
            :class:`~acalib.Container` with the cube slices, segmentated images and region of interest tables for each scale analyzed.
        """

        if isinstance(cube,NDData):
            if cube.wcs:
                wcs = cube.wcs
            else:
                wcs = None
            data = cube.data
        else:
            data = cube
            wcs = None


        c = acalib.Container()
        params = {"P":self.config["P"], "PRECISION":self.config["PRECISION"]}
        gms = GMS(params)


        spectra, slices = acalib.core.spectra_sketch(data, self.config["SAMPLES"], self.config["RANDOM_STATE"])

        pp_slices = []
        for slice in slices:
            pp_slice = vel_stacking(cube, slice)
            labeled_images = gms.run(pp_slice)

            if wcs is not None:
                freq_min = float(wcs.all_pix2world(0, 0, slice.start, 1)[2])
                freq_max = float(wcs.all_pix2world(0, 0, slice.stop, 1)[2])
            else:
                freq_min = None
                freq_max = None

            table = measure_shape(pp_slice, labeled_images, freq_min, freq_max)
            if len(table) > 0:
                c.tables.append(table)
                c.images.append(pp_slice)
                c.images.extend(labeled_images)

        if wcs:
            wcs = wcs.dropaxis(2)
            for i,im in enumerate(c.images):
                c.images[i] = NDDataRef(data=im, wcs = wcs)

        c.images.insert(0, cube)
        c.primary = c.images[0]
        return c
