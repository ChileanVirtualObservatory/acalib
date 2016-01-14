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
from mayavi import mlab

class BubbleClumps:

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
      self.par['SNRLIMIT']=1.0
      # Data Reduction Limit
      self.par['MAXBUB']=1000
      # Bubble Size = How many resolution FHWMs to consider
      self.par['BSIZE']=1.0
      # Cut level (how many RMSs we consider to force the compact support)
      self.par['CUTLEV']=1.0

   def _create_bubble(self):
      """This function creates a sub-cube containing a bubble, using the already computed values of 
         $\Delta_b$, $\Delta_s$, $\sigma_b$, and $\sigma_s$.
          
         $bub(p)=\exp(-\frac{(p - \Delta/2)^\top \Lambda (p - \Delta/2)}{2.0})$
         $p = [z y x]^\top$
         $\Delta = [\Delta_s \Delta_b \Delta_b]^\top$
         $\Lambda = \begin{array}{ccc} 1/\sigma_s^2 & 0 & 0 \\ 0 & 1/\sigma_s^2 & 0 \\ 0 & 0 & 1/\sigma_b^2  \end{array}$

         returns: a 3D numpy array containing the bubble (bub)
      """
      ds=self.ds
      db=self.db
      x=np.arange(-ds,ds+1)
      y=np.arange(-db,db+1)
      z=np.arange(-db,db+1)
      xyz=np.meshgrid(x,y,z,indexing='ij')
      ii=np.empty((3,len(x)*len(y)*len(z)))
      ii[0]=xyz[0].ravel()/self.ss
      ii[1]=xyz[1].ravel()/self.sb
      ii[2]=xyz[2].ravel()/self.sb
      V=(ii*ii)
      quad=V.sum(axis=0)
      res=np.exp(-quad/2.0)
      val=(res/res.max())
      bub=val.reshape(2*ds+1,2*db+1,2*db+1)
      return bub

   def _eighth_bubble(self):
      """This function creates an eighth of a sub-cube containing a bubble.
         It works similar to _create_bubble but it returns the values
         in value-features arrays. The eighth of a bubble is usefull because 
         the symetry of Gaussians (less values to process).

         returns: 
            - val: a 1D numpy array with the intensity values
            - idx: a 2D numpy array with the positions of these values in z, y, and x
      """
      ds=self.ds
      db=self.db
      x=np.arange(0,ds+1)
      y=np.arange(0,db+1)
      z=np.arange(0,db+1)
      xyz=np.meshgrid(x,y,z,indexing='ij')
      ii=np.empty((3,len(x)*len(y)*len(z)))
      idx=np.empty((3,len(x)*len(y)*len(z)))
      idx[0]=xyz[0].ravel()
      idx[1]=xyz[1].ravel()
      idx[2]=xyz[2].ravel()
      ii[0]=idx[0]/self.ss
      ii[1]=idx[1]/self.sb
      ii[2]=idx[2]/self.sb
      V=(ii*ii)
      quad=V.sum(axis=0)
      res=np.exp(-quad/2.0)
      val=(res/res.max())
      return val,idx

   def _update_min_energy(self,mat,ub,lb,delta):
     """Updates the minimum energies of self.energy from mat defaced by delta. 
        ub and lb bounds are provided to shrink the mat matrix when required (out of bounds, or partial update)
     """
     bord=np.array(self.energy.shape())
     # TODO: Bad usage of AData, because we want a reference of the data to modify it
     ene=self.energy.data
     # Numpyfy everithing
     ub=np.array(ub)
     lb=np.array(lb)
     delta=np.array(delta,dtype=int)
     # Create energy (e) and mat (m) indices 
     eub=ub + delta
     elb=lb + delta
     mub=ub-lb
     mlb=np.array([0,0,0])
     umask=eub > bord
     lmask=elb < 0
     mub[umask]-=eub[umask] - bord[umask]
     mlb[lmask]-=elb[lmask]
     eub[umask]=bord[umask]
     elb[lmask]=0
     # Obtain a reduced view of the matrices
     eview=ene[elb[0]:eub[0],elb[1]:eub[1],elb[2]:eub[2]]
     mview=mat[mlb[0]:mub[0],mlb[1]:mub[1],mlb[2]:mub[2]]
     # Select those that are lower in mat than in energy
     cmat=mview < eview
     # Update them in the energy matrix.
     eview[cmat]=mview[cmat]
     
   def _update_energies(self,lb,ub):
      """Update the energies, only from the lb to the ub points. 
      """
      #TODO: I do now know if get_slice actually states that we are making a copy...
      mcb=self.residual.get_slice(lb,ub)
      #Obtain the reference of the eighth of the bubble.
      vv=self.eival
      ff=self.eifeat.T
      # Iterates for every point in the eighth of the bubble
      # this starts from one because we do not want to repeat the position (0,0,0)
      for i in range(1,vv.size):
         mat=mcb/vv[i]
         d=ff[i]
         # update in the eight directions
         delta=np.array([1,1,1])*d
         self._update_min_energy(mat,ub,lb,delta)
         delta=np.array([1,1,-1])*d
         self._update_min_energy(mat,ub,lb,delta)
         delta=np.array([1,-1,1])*d
         self._update_min_energy(mat,ub,lb,delta)
         delta=np.array([1,-1,-1])*d
         self._update_min_energy(mat,ub,lb,delta)
         delta=np.array([-1,1,1])*d
         self._update_min_energy(mat,ub,lb,delta)
         delta=np.array([-1,1,-1])*d
         self._update_min_energy(mat,ub,lb,delta)
         delta=np.array([-1,-1,1])*d
         self._update_min_energy(mat,ub,lb,delta)
         delta=np.array([-1,-1,-1])*d
         self._update_min_energy(mat,ub,lb,delta)

  
   def fit(self,cube,verbose=False,use_meta=True):
      """Feed the algorithm with a Cube to bubblelize. This process generates:
         - A synthetic cube with the subtracted values (self.syn)
         - A residual cube with the original cube minus the synthetic (self.residual)
         - The positions of the fitted bubbles (self.positions)
         - The omplitued of the fitted bubbles (self.amplitudes)
         - The energy cube, which is aneroded version of the residual (self.energy) """


      # Set the RMS, or automatically find an estimate for it
      if not self.par.has_key('RMS'):
         rms=cube.estimate_rms()
         self.par['RMS']=rms
      rms=self.par['RMS']
      #plt.ion()
      # TODO: set parameters according to meta
      velres=self.par['VELORES']
      fwhmbeam=self.par['FWHMBEAM']
      bsize=self.par['BSIZE']
      snrlimit=self.par['SNRLIMIT']
      maxbub=self.par['MAXBUB']
      cutlev=self.par['CUTLEV']
      # Copy the supplied cube into a work cube which will hold the
      # residuals remaining after subtraction of the fitted Gaussians. 
      self.orig=cube
      self.residual=cube.copy()
      datamax=self.residual.max()
      datamax=datamax[0]
      self.syn=cube.empty_like()
      self.energy=cube.copy()
      #TODO: wrong usage of AData. I need to fill all non-nan values with a value (datamax)
      mas=np.isnan(self.energy.data)
      self.energy.data[np.logical_not(mas)]=datamax
      # Sigma values
      self.sb=bsize*fwhmbeam*self.FWHM_TO_SIGMA
      self.ss=bsize*velres*self.FWHM_TO_SIGMA
      self.S=np.array([[1.0/np.square(self.ss),0,0],[0,1.0/np.square(self.sb),0],[0,0,1/np.square(self.sb)]])
      # Deltas
      self.db=int(np.sqrt(-2*self.sb*self.sb*np.log(cutlev*rms/datamax)))
      self.ds=int(np.sqrt(-2*self.ss*self.ss*np.log(cutlev*rms/datamax)))
      if verbose:
         log.info("Datamax =="+str(datamax))
         log.info("RMS ="+str(rms))
         log.info("Computed Deltas ="+str((self.db,self.ds)))
      (self.eival,self.eifeat)=self._eighth_bubble()
      lb=(0,0,0)
      ub=self.residual.shape()
      self._update_energies(lb,ub)
      cb=self._create_bubble()
      delta=np.array([self.ds,self.db,self.db])
      iterate = True
      niter=0
      self.amplitudes=[]
      self.positions=[]
      while True:
         # Report the iteration number to the user if required.
         niter+=1
         if (niter > maxbub):
            if verbose:
               log.info("Criterion Met: maximum bubbles = "+str(maxbub))
            break
         y,xmax=self.energy.max()
         xmax=np.array(xmax)
         rem=y - rms
         if rem <= 0.0:
            if verbose:
               log.info("Criterion Met: Energy == 0 ")
            break
         self.amplitudes.append(rem)
         self.positions.append(xmax)
         if (y/rms - 1.0 < snrlimit):
            if verbose:
               log.info("Criterion Met: SNR="+str(y/rms-1.0)+"<"+str(snrlimit))
            break
         if verbose:
            log.info("Iteration: "+str(niter))
            log.info("Maximum energy E = "+str(y)+" at "+str(xmax))
            log.info("Remove E = "+str(rem)+" SNR="+str(y/rms - 1.0))
         ub=xmax + delta + 1
         lb=xmax - delta
         self.residual.add_flux(-rem*cb,lb,ub)
         self.syn.add_flux(rem*cb,lb,ub)
         self._update_energies(lb,ub)
      # Type hack..
      self.positions=np.array(self.positions)
      # TODO: put this in cube
      sh=cube.data.shape
      xi, yi, zi = np.mgrid[0:sh[0], 0:sh[1], 0:sh[2]]
      print xi.shape
      print cube.data.shape
      grid = mlab.pipeline.scalar_field(xi, yi, zi, cube.data)
      mmin = cube.data.min()
      mmax = cube.data.max()
      figure = mlab.figure('Orig')
      mlab.pipeline.volume(grid)
#, vmin=mmin, vmax=mmin + .*(mmax-mmin))
      mlab.axes()
      #mlab.show()
      figure = mlab.figure('Synthetic')
      grid = mlab.pipeline.scalar_field(xi, yi, zi, self.syn.data)
      mmin = self.syn.data.min()
      mmax = self.syn.data.max()
      mlab.pipeline.volume(grid)
#, vmin=mmin, vmax=mmin + .5*(mmax-mmin))
      mlab.axes()
      mlab.show()

      

   def linkage(self,verbose=False,show_dendogram=False):
      vect=self.positions
      print vect
      self.D=dist.pdist(vect,lambda x,y: np.exp(-(x-y).T*self.S*(x-y)))
      self.link=hier.linkage(D)
      if show_dendogram:
         hier.dendrogram(self.link)
         plt.show()

   def reasonable_cluster(self):
       two_delta=2*self.db*np.sqrt(3)
       cl=self.clustering(two_delta)
       plt.clf()
       fig = plt.figure(1, figsize=(4, 3))
       print " eps="+str(two_delta)+" k="+str(cl.labels_.max()+1)
       self.draw_cluster(fig,cl)
       plt.show()

   def selected_clusters(self,n_sols):
       v_ini=1.0
       self.solution=[]
       #self.scores=[]
       cl=self.clustering(1.0)
       old_labels=cl.labels_
       labels=old_labels
       while True:
          cl=self.clustering(v_ini)
          labels=cl.labels_
          #score=ami_score(labels,old_labels)
          #self.scores.append(score)
          lab_max=labels.max()
          lab_min=labels.min()
          print v_ini,lab_max,lab_min
          old_labels=labels
          v_ini+=1
          if lab_max > -1:
             self.solution.append(cl)
          else:
             continue
          if lab_min == lab_max:
              break
       #print self.scores
       #plt.clf()
       #plt.plot(self.scores)
       #plt.show()
       si=len(self.solution)
       affin=np.zeros((si,si))
       for i in range(si):
          for j in range(si):
             affin[i,j]=ami_score(self.solution[i].labels_,self.solution[j].labels_)
       affin+=1
       sclust=SpectralClustering(n_clusters=n_sols,affinity='precomputed')
       sclust.fit(affin)
       slabels=sclust.labels_
       sscore=[0]*n_sols
       sol_tmp=[None]*n_sols
       sol_eps=[0]*n_sols
       for i in range(si):
          k=slabels[i]
          c_score=affin[i,slabels==k].sum()
          if c_score > sscore[k]:
             sscore[k]=c_score
             sol_tmp[k]=self.solution[i]
             sol_eps[k]=i
       self.solution=sol_tmp
       plt.clf()
       f_num=0
       for cl in self.solution:
           f_num+=1
           fig = plt.figure(f_num, figsize=(4, 3))
           print "sol "+str(f_num)+" eps="+str(sol_eps[f_num-1])+" k="+str(cl.labels_.max()+1)
           self.draw_cluster(fig,cl)
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

   def agglo_clust(self,val,criterion='inconsistent'):
       lab=hier.fcluster(self.link,val,criterion=criterion)
       return namedtuple('labels_',lab,'linkage_',self.link,'pdist_',self.dist)

   def draw_cluster(self,fig,clust):
       
       ax = Axes3D(fig, rect=[0, 0, .95, 1], elev=48, azim=134)

       plt.cla()
       amp=np.array(self.amplitudes)
       amp=300.0*amp/amp.max()
       labels = clust.labels_
       lab_select=labels[labels>-1]
       pos_select=self.positions[labels>-1,:]
       sur_select=amp[labels>-1]
       pos_noselect=self.positions[labels==-1,:]
       sur_noselect=amp[labels==-1]
       ax.scatter(pos_noselect[:, 2], pos_noselect[:, 1], pos_noselect[:, 0], c='black',alpha=0.5,s=sur_noselect)
       ax.scatter(pos_select[:, 2], pos_select[:, 1], pos_select[:, 0], c=lab_select.astype(np.float),alpha=0.5,s=sur_select)

       ax.w_xaxis.set_ticklabels([])
       ax.w_yaxis.set_ticklabels([])
       ax.w_zaxis.set_ticklabels([])
       ax.set_xlabel('Declination')
       ax.set_ylabel('Right Ascension')
       ax.set_zlabel('Velocity')


   def test_clustering(self):

       plt.clf()
       cl=self.clustering(30,'dbscan')
       fig = plt.figure(1, figsize=(4, 3))
       plt.title('dbscan')
       self.draw_cluster(fig,cl)

       cl=self.clustering(0.7,'affinity_propagation')
       fig = plt.figure(2, figsize=(4, 3))
       plt.title('affinity_propagation')
       self.draw_cluster(fig,cl)

       cl=self.clustering(5,'kmeans')
       fig = plt.figure(3, figsize=(4, 3))
       plt.title('kmeans')
       self.draw_cluster(fig,cl)

       cl=self.clustering(5,'spectral')
       fig = plt.figure(4, figsize=(4, 3))
       plt.title('spectral')
       self.draw_cluster(fig,cl)

       #self.linkage()
       #cl=self.clustering(0.8,'agglomerative')
       #fig = plt.figure(5, figsize=(4, 3))
       #plt.title('agglomerative')

       #self.draw_cluster(fig,cl)
       plt.show()
 


      
#   def recursive_candidate(self,l_labels,r_labels,val,l_val,r_val,ami_th=0.6,granularity=1.0,method='dbscan'):
#       if r_val - l_val < granularity:
#          return
#       clust=self.clustering(val,method)
#       labels=clust.labels_
#       print labels
#       print r_labels
#       print l_labels
#       l_ami=ami_score(labels,l_labels)
#       r_ami=ami_score(labels,r_labels)
#       #print l_labels
#       #print r_labels
#       print "val = "+str(val)
#       #print "l_val = "+str(l_val)
#       #print "r_val = "+str(r_val)
#       #print "l_ami = "+str(l_ami)
#       #print "r_ami = "+str(r_ami)
#       flag=0;
#       if (l_ami < ami_th):
#          flag=1
#          self.recursive_candidate(l_labels,labels,(l_val + val)/2.0,l_val,val,ami_th,granularity,method)
#       if (r_ami < ami_th):
#          flag=1
#          self.recursive_candidate(labels,r_labels,(r_val + val)/2.0,val,r_val,ami_th,granularity,method)
#       if flag!=0:
#          self.solution.append(clust)
#          print "candidates = "+str(len(self.solution))
#
#   def cluster_candidates(self,ami_th=0.6,method='dbscan'):
#       pos=self.positions
#       if method=='dbscan':
#          v_ini=2.0
#          v_end=self.positions.max()/2.0
#          granularity=1.0
#       elif method=='kmeans':
#          v_ini=1
#          v_end=self.positions.shape[0]/2.0
#          granularity=1.0
#       elif method=='agglomerative':
#          v_ini=0.0
#          v_end=10.0
#          granularity=0.01
#          clust=self.agglo_clust(val); # inconsistency
#       elif method=='affinity_propagation':
#          v_ini=0.5
#          v_end=0.9
#          granularity=0.001
#       elif method=='spectral':
#          v_ini=1
#          v_end=self.positions.shape[0]/2.0
#          granularity=1.0
#       log.info("v_ini="+str(v_ini)+" v_end="+str(v_end))
#       self.solution=[]
#       cl=self.clustering(v_ini,method)
#       self.solution.append(cl)
#       l_labels=cl.labels_
#       cl=self.clustering(v_end,method)
#       self.solution.append(cl)
#       r_labels=cl.labels_
#       self.recursive_candidate(l_labels,r_labels,(v_ini + v_end)/2.0,v_ini,v_end,ami_th=ami_th,granularity=granularity,method=method)
#       f_num=0
#       plt.clf()
#       for cl in self.solution:
#           f_num+=1
#           fig = plt.figure(f_num, figsize=(4, 3))
#           self.draw_cluster(fig,cl)
#       plt.show()


#   def clusterize(self,verbose=False):
#      # Binary Search
#      
#      pos=np.array(self.positions)
#      for j in range(20):
#         db = DBSCAN(eps=1.0+10.0*j,).fit(pos)
#         vect=db.labels_
#         k=int(np.max(db.labels_))
#         print vect
#         print k
#         colors = iter(cm.rainbow(np.linspace(0, 1, k+1)))
#         plt.imshow(self.syn.get_stacked(axis=0),cmap='Greys')
#         v=pos[vect==-1]
#         plt.scatter(v[:,2],v[:,1], color='black',marker='+',alpha=0.5)
#         for i in range(k+1):
#            v=pos[vect==i]
#            plt.scatter(v[:,2],v[:,1], color=next(colors),marker='o',alpha=0.5)
#         plt.show()
#
#
#   def clusterize_hier(self,verbose=False):
#      """Under development """
#      # Heriarchical Clustering
#      # Compute the condensated eucledian distance matrix
#      vect=np.array(self.positions)
#      M=dist.pdist(vect,lambda x,y: np.exp(-(x-y).T*self.S*(x-y)))
#      Z=hier.linkage(M)
#      hier.dendrogram(Z)
#      plt.show()
#      T=hier.fcluster(Z,0.0)
#      lim=float(max(T))/2.0
#      print "lim",lim
#      step=0.01
#      val=0.0
#      while True:
#         #print "val",val
#         val+=step
#         T=hier.fcluster(Z,val)
#         k=max(T)
#         #print "k,lim",k,lim
#         if k < lim:
#            break
#
#      colors = iter(cm.rainbow(np.linspace(0, 1, k+1)))
#      plt.imshow(self.syn.get_stacked(axis=0),cmap='Greys')
#      for i in range(k+1):
#         v=vect[T==i]
#         print v
#         plt.scatter(v[:,2],v[:,1], color=next(colors),marker='o',alpha=0.5)
#      plt.show()
#
#     
#
#
#
#
