import numpy as np
from astropy import log
from sklearn.cluster import DBSCAN,KMeans
import matplotlib.pyplot as plt
import acalib

_clustering_method_dict = {'DBSCAN': lambda param, rep: DBSCAN(eps=param).fit(rep),
               'KMEANS': lambda param, rep: KMeans(n_clusters=param).fit(rep)}

from graphviz import Digraph
import colorsys

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
            
    def synthesize(self,clumps):
        newSyn = np.zeros(self.synth.data.shape)
        if not isinstance(clumps,list):
            clumps=[clumps]
        for i in clumps:
            newSyn = acalib.core.synthesize_bubbles(newSyn, self.clumps[i], self.kernel, self.noise, self.delta)
        return acalib.upi.Data(newSyn * self.scale - self.shift, wcs=self.synth.wcs, unit=self.synth.unit, meta=self.synth.meta,mask=self.synth.mask)

    def _hide_all(self,parent,node):
        self.hide(node)
        if parent[node] is None:
            return
        for key in parent[node]:
            self._hide_all(parent[node],key)

    def contract(self,node):
        parent=self.find_parent(node, self.tree)
        self.unhide(node)
        if parent[node] is None:
            return
        for key in parent[node]:
            self._hide_all(parent[node],key)

    def expand(self,node):
        parent=self.find_parent(node, self.tree)
        if parent[node] is None:
            log.error("Node "+str(node)+" have no children")
            return
        self.hide(node)
        for key in parent[node]:
            self.unhide(key)

    def ungroup(self,node):
        self.expand(node)
        parent=self.find_parent(node, self.tree)
        stree=parent[node]
        for i in stree:
            parent[i]=stree[i]
        del parent[node]
        del self.expand_method[node]
        del self.display[node]
        del self.enabled[node]
        del self.clumps[node]

    def group(self,clumps):
        lst=[]
        for i in clumps:
            lst.append(self.find_parent(i,self.tree))
        if not lst[1:] == lst[:-1]:
            log.error("Clumps must have the same parent node to group them")
            return
        parent=lst[0]
        parent[self.nextnode] = dict()
        new_node=parent[self.nextnode]
        new_rep=[]
        for i in clumps:
            new_rep.extend(self.clumps[i])
            new_node[i]=parent[i]
            del parent[i]
            self._hide_all(new_node,i)
        self.expand_method[self.nextnode] = ("Manual",len(lst))
        self.display[self.nextnode] = True
        self.enabled[self.nextnode] = True
        self.clumps[self.nextnode] = new_rep
        self.nextnode += 1

    def add_to_group(self,node,clumps):
        parent_node=self.find_parent(node, self.tree)
        lst=[]
        for i in clumps:
            if parent_node != self.find_parent(i,self.tree):
                log.error("Clumps must have the same parent node to group them")
                return
        new_rep  = self.clumps[node]
        new_node = parent_node[node]
        for i in clumps:
            new_rep.extend(self.clumps[i])
            new_node[i]=parent_node[i]
            del parent_node[i]
            self._hide_all(new_node,i)

    def clusterize(self,node,method='KMEANS',param=2):
        parent=self.find_parent(node, self.tree)
        if parent is None:
            log.error("No such node "+str(node)+" in the tree")
            return
        if parent[node] is not None:
            log.error("Node "+str(node)+" already expanded, please remove the expansion first (tree.contract)")
            return
        if self.enabled[node] is False:
            log.error("Node " + str(node) + " is disabled, please enable it first (tree.enable)")
            return
        mfunc = _clustering_method_dict[method]
        myrep = self.clumps[node]
        clust = mfunc(param,myrep)
        imax = max(clust.labels_)
        parent[node]=dict()
        self.display[node] = False
        self.expand_method[node] = (method,param)
        subtree=parent[node]
        for i in range(imax+1):
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

    def _unclusterize_subtree(self, tree):
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

    def unclusterize(self,node):
        parent = self.find_parent(node, self.tree)
        if parent[node] is None:
            log.error("Node " + str(node) + " not expanded, cannot contract")
            return
        del self.expand_method[node]
        self._unclusterize_subtree(parent[node])
        self.display[node] = True

    def find_parent(self, node, tree):
        for elm in tree:
            if elm == node:
                return tree
            elif tree[elm] is not None:
                par = self.find_parent(node,tree[elm])
                if par is not None:
                    return par
        return None

    def show_tree(self,node=0):
        parent=self.find_parent(node,self.tree)
        n_colors = self.visible_clumps(parent,node)
        colors = plt.cm.rainbow(np.linspace(0, 1, n_colors))
        return self._create_digraph(node,colors)

    def show(self,node=0):
        self.show_clumps(node=node)
        return self.show_tree(node=node)

    def show_clumps(self,node=0):
        acalib.show_clumps(self,node=node)

    def visible_clumps(self,subtree,node):
        n=0
        if self.enabled[node]:
            if self.display[node]:
                n=1
        return self._visible_clumps(subtree[node]) + n

    def _visible_clumps(self,subtree):
        n = 0
        if subtree is not None:
            for elm in subtree:
                if self.enabled[elm]:
                    if self.display[elm]:
                        n += 1
                    n+=self._visible_clumps(subtree[elm])
        return n

    def _create_digraph(self, node,colors):
        dot = Digraph(comment='Tree')
        self.t_n=0
        parent=self.find_parent(node,self.tree)
        def recu(parent, tree):
            for (key, val) in tree.items():
                if self.enabled[key]:
                    if self.display[key]:
                        hsv = colorsys.rgb_to_hsv(colors[self.t_n][0], colors[self.t_n][1], colors[self.t_n][2])
                        dot.node(str(key),str(key),{'color': "%f, %f, %f" % (hsv[0], hsv[1], hsv[2]), 'style': 'filled'})
                        self.t_n+=1
                    dot.node(str(key), str(key))
                    if isinstance(val, dict):
                        recu(key, val)
                    if parent != key:
                        dot.edge(str(parent), str(key))
        if self.enabled[node]:
            if self.display[node]:
                hsv = colorsys.rgb_to_hsv(colors[self.t_n][0], colors[self.t_n][1], colors[self.t_n][2])
                dot.node(str(node),str(node),{'color': "%f, %f, %f" % (hsv[0], hsv[1], hsv[2]), 'style': 'filled'})
                self.t_n+=1
        if isinstance(parent[node], dict):
            recu(node, parent[node])
        return dot


