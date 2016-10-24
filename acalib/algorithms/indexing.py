import acalib
from .algorithm import Algorithm


class Indexing(Algorithm):

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
            pp_slice = acalib.core.vel_stacking(data.data, slice)
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
