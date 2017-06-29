import acalib
from .algorithm import Algorithm
from astrodendro import Dendrogram

class Astrodendro(Algorithm):

    def default_params(self):
       self.config.update({"min_value":"min", "min_delta":0, "min_npix":0,"is_independent":None, "verbose":False, "neighbours":None,"wcs":None})

    def run(self,data):
        self.dendrogram = Dendrogram.compute(data, **self.config)
        return self.dendrogram

    def data(self):
        return self.dendrogram.data

    def trunk(self):
        return self.dendrogram.trunk

    def neighbours(self, idx):
        return self.dendrogram.neighbours(idx)

    def leaves(self):
        return self.dendrogram.leaves()

    def structure_at(self, indices):
        return self.dendrogram.structure_at(indices)

    def all_structures(self):
        return self.dendrogram.all_structures()

    def plotter(self):
        return self.dendrogram.plotter()

    def show(self):
        self.dendrogram.viewer().show()
        return

    def index_map(self):
        return self.dendrogram.index_map