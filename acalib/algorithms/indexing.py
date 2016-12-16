import acalib
from .algorithm import Algorithm
from .gms import GMS


class Indexing(Algorithm):

    """
    Perform an unsupervised region of interest detection and extract shape features.

    Parameters
    ----------
    params: dict (default = None)
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

    def run(self, data):
        """
            Run the indexing algorithm on a given data cube.

            Parameters
            ----------            
            data : (M,N,Z) numpy.ndarray or astropy.nddata.NDData
                Astronomical data cube.

            Returns
            -------
            result: :class:`~acalib.Container` with the cube slices, segmentated images and region of interest tables for each scale analyzed.
        """

        if data.wcs:
            wcs = data.wcs
        else:
            wcs = None

        c = acalib.Container()
        params = {"P":self.config["P"], "PRECISION":self.config["PRECISION"]}
        gms = GMS(params)


        spectra, slices = acalib.core.spectra_sketch(data.data, self.config["SAMPLES"], self.config["RANDOM_STATE"])

        pp_slices = []
        for slice in slices:

            pp_slice = acalib.core.vel_stacking(data, slice)
            labeled_images = gms.run(pp_slice)

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
