import numpy as np
from astropy import log
from sklearn.cluster import DBSCAN,KMeans

import acalib

_clustering_method_dict = {'DBSCAN': lambda (param, rep): DBSCAN(eps=param).fit(rep),
               'KMEANS': lambda (param, rep): KMeans(n_clusters=param).fit(rep)}


class HRTree():


    def __init__(self,table,synth):
        self.nextnode = 0

        if table.meta['KERNEL'] != 'METABUBBLE':
            log.error("Only METABUBBLE representation supported for now")
            return
        self.rep, self.delta, self.noise, self.kernel = acalib.HRep.toTuple(table)
        self.synth = synth
        self.scale = table.meta['SCALE']
        self.shift = table.meta['SHIFT']
        self.clumps = dict()
        self.display = dict()
        self.enabled = dict()
        self.tree = dict()
        self.expand_method = dict()

        self.clumps[self.nextnode] = self.rep
        self.tree[self.nextnode] = None
        self.display[self.nextnode] = True
        self.enabled[self.nextnode] = True
        self.nextnode += 1


    def expand(self,node,method='KMEANS',param=2):
        parent=self.find_dict(node, self.tree)
        if parent == None:
            log.error("No such node "+str(node)+" in the tree")
            return
        if parent[node] != None:
            log.error("Node "+str(node)+" already expanded, please remove the expansion first (tree.contract)")
            return
        if self.enabled[node] == False:
            log.error("Node " + str(node) + " is disabled, please enable it first (tree.enable)")
        mfunc = _clustering_method_dict[method]
        myrep = self.clumps[node]
        clust=mfunc(param,myrep)
        ###expand

    def disable(self,node):
        pass

    def enable(self,node):
        pass

    def hide(self,node):
        pass

    def reveañ(self,node):
        pass

    def contract(self,node):
        pass
