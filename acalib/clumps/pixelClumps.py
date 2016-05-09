import numpy as np
from collections import namedtuple
import copy
import matplotlib.pyplot as plt
import sys
from astropy import log
import scipy.spatial.distance as dist
import scipy.cluster.hierarchy as hier
from sklearn.cluster import DBSCAN
from sklearn.cluster import KMeans
from sklearn.cluster import AffinityPropagation
from sklearn.cluster import SpectralClustering
from sklearn.metrics.cluster import adjusted_mutual_info_score as ami_score
from mpl_toolkits.mplot3d import Axes3D
from sklearn import metrics
import matplotlib.cm as cm
from scipy import ndimage
from scipy.spatial.distance import pdist
import math
from acalib.io import graph
from acalib import flux
from numpy.linalg import inv
from numpy.linalg import det
from numpy.linalg import norm
from compiler.ast import flatten

class PixelClumps:

   def __init__(self):
      self.FWHM_TO_SIGMA = 1. / (8 * np.log(2))**0.5
      # Set the very 
      self.defaultParams()
    
   def defaultParams(self):
      self.par=dict()
      # Spectral Resoluion in pixels (smoothing function)
      self.par['VELORES']=2.0
      # Beam resoluion in pixels (smoothing function)
      self.par['FWHMBEAM']=2.0
      # How many RMS is considered noise
      self.par['SNRLIMIT']=2.0
      # Agglomerative Criterion
      self.par['AGGLOCRIT']='distance'
      # Agglomerative Parameter Limit (stoppong criterion)
      self.par['AGGLOLIMIT']=4.0

   #def add_element(self,new_amp,new_mu,new_Sigma)
   #   n=self.amps.shape[0]
      
   def join_gauss(self,a1,a2,x1,x2,P1,P2):
      m1=a1*np.sqrt(8.0*np.pi*np.pi*np.pi/det(P1))
      m2=a2*np.sqrt(8.0*np.pi*np.pi*np.pi/det(P2))
      m=m1 + m2
      my_mu=(m1*x1 +m2*x2)/m
      S12=np.outer((x1 - x2),(x1 - x2))
      my_Sigma=(m1*inv(P1) + m2*inv(P2) + m1*m2*S12/m)/m
      my_P=inv(my_Sigma)
      my_a=m/(np.sqrt(8.0*np.pi*np.pi*np.pi/det(my_P)))
      
      return my_a,my_mu,my_P
 
   def KL_div(self,pixel_tree,new_mu,new_Prec):
         # Evaluate
         elems=flatten(pixel_tree)
         e_amps=self.orig_amps[elems]
         e_amps=e_amps/e_amps.sum()
         feat=self.orig_mu[elems]
         gauss=flux.create_gauss(new_mu,new_Prec,feat.T,1.0)
         gauss=gauss/gauss.sum()
         KL=e_amps.dot(np.log2(e_amps/gauss))
         return(KL)

   def ISD(self,i,j):
         pass 

   
   def radial_overlap(self,i,j):
         rms=self.par['RMS']
         v=self.mu[i] - self.mu[j]
         u=v/norm(v)
         q1=np.dot(u,np.dot(self.Prec[i],u))
         q2=np.dot(u,np.dot(self.Prec[j],u))
         r1=np.sqrt(2*np.log(self.amps[i]/rms)/q1)
         r2=np.sqrt(2*np.log(self.amps[j]/rms)/q2)
         return r1 + r2 - norm(v)

   #def KL_vector(self,elem):
   #      n=self.amps.shape[0]
   #      vect=np.zeros(n)
   #      for j in range(n):
   #          (new_amp,new_mu,new_Prec)=self.join_gauss(self.amps[elem],self.amps[j],self.mu[elem],self.mu[j],self.Prec[elem],self.Prec[j])
   #          vect[j]=self.KL_div((self.pixel_tree[elem],self.pixel_tree[j]),new_mu,new_Prec)
   #      return vect

   def compare_to_all(self,i):
         if self.par['AGGLOCRIT']=='distance':
            darr=dist.cdist([self.mu[i]],self.mu)
            darr=darr[0]
            return darr
         if self.pas['AGGLOCRIT']=='radial_overlap':
            darr=np.zeros(self.amps.shape)
            for j in range(self.amps.shape):
               darr[j]=self.radial_overlap(i,j)
         if self.pas['AGGLOCRIT']=='kl_divergence':
               pass
         return darr

   def fit(self,cube,verbose=False,use_meta=True):
      # Set the RMS, or automatically find an estimate for it
      if not self.par.has_key('RMS'):
         rms=cube.estimate_rms()
         self.par['RMS']=rms
      rms=self.par['RMS']
      self.orig=cube
      rms=cube.estimate_rms()
      threshold=self.par['SNRLIMIT']*rms
      ff=np.where(cube.data>threshold)
      self.amps=cube.data[ff].filled()
      self.mu=np.transpose(ff).astype(float)
      self.orig_amps=self.amps
      self.orig_mu=self.mu
      I=np.identity(3)
      self.Prec=[]
      for a in self.amps:
          val=2*np.log(a/rms)
          self.Prec.append(val*I)
      self.Prec=np.array(self.Prec)
      n=self.amps.shape[0]
      self.pixel_tree=[]
      minimals=np.zeros(self.amps.shape)
      mindices=np.zeros(self.amps.shape)
      for i in range(0,n):
         self.pixel_tree.append((i))
         # Obtain performance vector
         darr=self.compare_to_all(i)
         darr[i]=np.Inf
         j=np.argmin(darr)
         mindices[i]=j
         minimals[i]=darr[j]
      while minimals.min() < self.par['AGGLOLIMIT'] and n > 2:
         print "iter",n
         print "min",minimals.min()
         i=int(np.argmin(minimals))
         j=int(mindices[i])
         (new_amp,new_mu,new_Prec)=self.join_gauss(self.amps[i],self.amps[j],self.mu[i],self.mu[j],self.Prec[i],self.Prec[j])
         i_pixels=self.pixel_tree[i]
         j_pixels=self.pixel_tree[j]
         if i>j:
             bak=i
             i=j
             j=bak
         # Remove elements
         self.amps=np.delete(self.amps,[i,j])
         self.mu=np.delete(self.mu,[i,j],axis=0)
         self.Prec=np.delete(self.Prec,[i,j],axis=0)
         mindices=np.delete(mindices,[i,j])
         minimals=np.delete(minimals,[i,j])
         self.pixel_tree.remove(i_pixels)
         self.pixel_tree.remove(j_pixels)
         # Add elements
         #pos=np.searchsorted(self.amps,new_amp)
         #pos=0
         #self.amps=np.insert(self.amps,pos,new_amp)
         #self.mu=np.insert(self.mu,pos,new_mu,axis=0)
         #self.Prec=np.insert(self.Prec,pos,new_Prec,axis=0)
         #mindices=np.insert(mindices,pos,-1) #Deleyed assignment for convinience
         #minimals=np.insert(minimals,pos,-1) #Deleyed assignment for convinience
         #self.pixel_tree.insert(pos,(i_pixels,j_pixels))

         # Select items to update
         recalc=np.hstack((np.where(mindices==j),np.where(mindices==i)))
         recalc=np.unique(recalc)
         # Correct all the other indices
         mask=np.where(mindices>j)
         mindices[mask]=mindices[mask]-1
         mask=np.where(mindices>i)
         mindices[mask]=mindices[mask]-1
         #mask=np.where(mindices>=pos)
         #mindices[mask]=mindices[mask]+1
         # Update Items
         for k in recalc:
            darr=self.compare_to_all(k)
            darr[k]=np.Inf
            l=np.argmin(darr)
            mindices[k]=l
            minimals[k]=darr[l]
         self.amps=np.append(self.amps,new_amp)
         self.mu=np.append(self.mu,[new_mu],axis=0)
         self.Prec=np.append(self.Prec,[new_Prec],axis=0)
         darr=self.compare_to_all(new_mu)
         darr[-1]=np.Inf
         l=np.argmin(darr)
         mindices=np.append(mindices,l)
         minimals=np.append(minimals,darr[l])
         self.pixel_tree.append((i_pixels,j_pixels))
         mask=np.where(darr < minimals)
         minimals[mask]=darr[mask]
         n=n-1
      #print self.amps
      #print self.mu
      #print self.Prec
      i=1
      labels=np.empty_like(self.orig_amps)
      for pix in self.pixel_tree:
          if type(pix) is int:
             elems=[pix]
          else:
             elems=flatten(pix)
          print "cluster",i,"=",len(elems)
          if len(elems) < 10:
             print "ignored!"
             labels[elems]=0
          else:
             labels[elems]=i
             i+=1
      self.draw_cluster(labels)
   
   def draw_cluster(self,labels):
       fig=plt.figure()
       ax = Axes3D(fig, rect=[0, 0, .95, 1], elev=48, azim=134)
       plt.cla()
       #amp=np.array(self.orig_amps)
       #amp=300.0*amp/amp.max()
       ax.scatter(self.orig_mu[:, 2], self.orig_mu[:, 1], self.orig_mu[:, 0], c=labels.astype(np.float),alpha=0.5)#,s=amp)
       ax.w_xaxis.set_ticklabels([])
       ax.w_yaxis.set_ticklabels([])
       ax.w_zaxis.set_ticklabels([])
       ax.set_xlabel('Declination')
       ax.set_ylabel('Right Ascension')
       ax.set_zlabel('Velocity')
       plt.show()



   def clustering(self,val,method='dbscan'):
       pos=self.positions
       if method=='dbscan':
          clust=DBSCAN(eps=val).fit(pos) # epsilon
       elif method=='kmeans':
          clust=KMeans(n_clusters=int(val)).fit(pos) # k-clusters
       elif method=='agglomerative':
          clust=self.agglo_clust(val); # inconsistency
       elif method=='affinity_propagation':
          clust=AffinityPropagation(damping=val).fit(pos)  # damping
       elif method=='spectral': 
          clust=SpectralClustering(n_clusters=int(val)).fit(pos) # n-clusters
       else:
          log.warning("The clustering algorithm is not supported yet.")
       return clust
