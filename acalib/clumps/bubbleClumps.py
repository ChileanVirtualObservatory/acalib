import numpy as np

import copy
import matplotlib.pyplot as plt
import sys
from astropy import log
import scipy.spatial.distance as dist
import scipy.cluster.hierarchy as hier
from sklearn.cluster import DBSCAN
from sklearn import metrics
import matplotlib.cm as cm


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
      self.S=2*np.array([[1.0/np.square(self.ss),0,0],[0,1.0/np.square(self.sb),0],[0,0,1/np.square(self.sb)]])
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

   def clusterize(self,verbose=False):
      pos=np.array(self.positions)
      for j in range(20):
         db = DBSCAN(eps=1.0+2.0*j,).fit(pos)
         vect=db.labels_
         k=int(np.max(db.labels_))
         print vect
         print k
         colors = iter(cm.rainbow(np.linspace(0, 1, k+1)))
         plt.imshow(self.syn.get_stacked(axis=0),cmap='Greys')
         v=pos[vect==-1]
         plt.scatter(v[:,2],v[:,1], color='black',marker='o',alpha=0.5)
         for i in range(k+1):
            v=pos[vect==i-1]
            plt.scatter(v[:,2],v[:,1], color=next(colors),marker='o',alpha=0.5)
         plt.show()
  
 
   def clusterize_hier(self,verbose=False):
      """Under development """
      # Heriarchical Clustering
      # Compute the condensated eucledian distance matrix
      vect=np.array(self.positions)
      M=dist.pdist(vect,lambda x,y: np.exp(-(x-y).T*self.S*(x-y)))
      Z=hier.linkage(M)
      hier.dendrogram(Z)
      plt.show()
      T=hier.fcluster(Z,0.0)
      lim=float(max(T))/2.0
      print "lim",lim
      step=0.01
      val=0.0
      while True:
         #print "val",val
         val+=step
         T=hier.fcluster(Z,val)
         k=max(T)
         #print "k,lim",k,lim
         if k < lim:
            break

      colors = iter(cm.rainbow(np.linspace(0, 1, k+1)))
      plt.imshow(self.syn.get_stacked(axis=0),cmap='Greys')
      for i in range(k+1):
         v=vect[T==i]
         print v
         plt.scatter(v[:,2],v[:,1], color=next(colors),marker='o',alpha=0.5)
      plt.show()

     




