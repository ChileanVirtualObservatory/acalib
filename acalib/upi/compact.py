import numpy as np
from astropy import log
from sklearn.cluster import DBSCAN,KMeans

import acalib

_clustering_method_dict = {'DBSCAN': lambda (param, rep): DBSCAN(eps=param).fit(rep),
               'KMEANS': lambda (param, rep): KMeans(n_clusters=param).fit(rep)}

from graphviz import Digraph

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
        parent=self.find_parent(node, self.tree)
        if parent is None:
            log.error("No such node "+str(node)+" in the tree")
            return
        if parent[node] is not None:
            log.error("Node "+str(node)+" already expanded, please remove the expansion first (tree.contract)")
            return
        if self.enabled[node] is False:
            log.error("Node " + str(node) + " is disabled, please enable it first (tree.enable)")
        mfunc = _clustering_method_dict[method]
        myrep = self.clumps[node]
        clust = mfunc(param,myrep)
        imax = max(clust.labels_)
        parent[node]=dict()
        self.display[node] = False
        self.expand_method[node] = (method,param)
        subtree=parent[node]
        for i in range(imax):
            new_rep = myrep[clust.labels_ == i]
            self.clumps[self.nextnode] = new_rep
            subtree[self.nextnode] = None
            self.display[self.nextnode] = True
            self.enabled[self.nextnode] = True
            self.nextnode += 1

    def disable(self,node):
        self.enabled[node] = False

    def enable(self,node):
        self.enabled[node] = True

    def hide(self,node):
        self.display[node] = False

    def unhide(self,node):
        self.display[node] = True

    def _contract_subtree(self, tree):
        if tree is None:
            return
        for elm in tree:
            if tree[elm] is not None:
                del self.expand_method[elm]
            self._contract_subtree(tree[elm])
            del self.display[elm]
            del self.enabled[elm]
            del self.clumps[elm]
            del tree[elm]

    def contract(self,node):
        parent = self.find_parent(node, self.tree)
        if parent[node] is None:
            log.error("Node " + str(node) + " not expanded, cannot contract")
            return
        del self.expand_method[node]
        self._contract_subtree(parent[node])
        self.display[node] = True

    def find_parent(self, node, tree):
        for elm in tree:
            if elm == tree:
                return tree
            elif tree[elm] is not None:
                par = self.find_parent(node,tree[elm])
                if par is not None:
                    return par
        return None

    def show(self):
        n_colors = self._reachable_clumps()
        # Create colors
        graphical_tree = self._create_digraph(n_colors)
        three_view_clumps = self._create_three_view(n_colors)
