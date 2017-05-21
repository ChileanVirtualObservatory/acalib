import acalib
from .algorithm import Algorithm
from .gms import GMS

import os
import distributed
import dask.bag as db
from astropy import log

class Indexing(Algorithm):
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

    def run(self, data):
        """
            Run the indexing algorithm on a given data cube.

            Parameters
            ----------
            data : (M,N,Z) numpy.ndarray or astropy.nddata.NDData
                Astronomical data cube.

            Returns
            -------
            :class:`~acalib.Container` with the cube slices, segmentated images and region of interest tables for each scale analyzed.
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

class IndexingDask(Algorithm):
    #Create master scheduler => dask-scheduler
    #Add slave to scheduler => dask-worker --nprocs [N Processes] [SCHEDULER_ADDR]
    def default_params(self):
        if 'P' not in self.config:
            self.config['P'] = 0.05
        if 'PRECISION' not in self.config:
            self.config['PRECISION'] = 0.02
        if 'RANDOM_STATE' not in self.config:
            self.config['RANDOM_STATE'] = None
        if 'SAMPLES' not in self.config:
            self.config['SAMPLES'] = 1000
        if 'N_PARTITIONS' not in self.config:
            self.config['N_PARTITIONS'] = None
        if 'PARTITION_SIZE' not in self.config:
            self.config['PARTITION_SIZE'] = None
        if 'SCHEDULER_ADDR' not in self.config:
            self.config['SCHEDULER_ADDR'] = '127.0.0.1:8786'

    def computeIndexing(self, data):
        gmsParams = {'P': self.config['P'], 'PRECISION': self.config['PRECISION']}
        gms = GMS(gmsParams)
        spectra, slices = acalib.core.spectra_sketch(data.data, self.config["SAMPLES"], self.config["RANDOM_STATE"])
        result = []
        for slice in slices:
            slice_stacked = acalib.core.vel_stacking(data, slice)
            labeled_images = gms.run(slice_stacked)
            freq_min = None
            freq_max = None
            if data.wcs:
                freq_min = float(data.wcs.all_pix2world(0, 0, slice.start, 1)[2])
                freq_max = float(data.wcs.all_pix2world(0, 0, slice.stop, 1)[2])
            table = acalib.core.measure_shape(slice_stacked, labeled_images, freq_min, freq_max)
            if len(table) > 0:
                result.append(table)
        return result

    def checkAbsoluteLocalFilePaths(self, files):
        for f in files:
            if not os.path.isabs(f):
                log.error('FITS file path should be absolute when running in local-filesystem mode')
                raise ValueError('FITS file path should be absolute when running in local-filesystem mode')

    def run(self, files):
        self.checkAbsoluteLocalFilePaths(files)
        log.info('Connecting to dask-scheduler at ['+self.config['SCHEDULER_ADDR']+']')
        client = distributed.Client(self.config['SCHEDULER_ADDR'])
        indexing = lambda x: self.computeIndexing(x)
        indexing.__name__ = 'computeIndexing'
        load = lambda x: acalib.io.loadFITS_PrimmaryOnly(x)
        load.__name__ = 'loadData'
        denoise = lambda x: acalib.denoise(x, threshold=acalib.noise_level(x))
        denoise.__name__ = 'denoise'
        cores = sum(client.ncores().values())
        log.info('Computing "Indexing" on '+str(len(files))+' elements with '+str(cores)+' cores')
        data = db.from_sequence(files, self.config['PARTITION_SIZE'], self.config['N_PARTITIONS'])
        data = data.map(load).map(denoise)
        results = data.map(indexing).compute()
        log.info('Gathering results')
        results = client.gather(results)
        log.info('Removing dask-client')
        client.shutdown()
        return results
