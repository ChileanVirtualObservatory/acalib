import numpy as np
from collections import deque
from scipy.optimize import fmin_bfgs,check_grad,approx_fprime
#from scipy.optimize.linesearch import (line_search_BFGS, line_search_wolfe1, line_search_wolfe2, line_search_wolfe2 as line_search)
#from optimize import fmin_bfgs

import copy
import matplotlib.pyplot as plt
import sys
from astropy import log


class PerfBubbleClumps:

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
      self.par['NRMS']=3.0
      # How meny FHWMs to consider
      self.par['BSIZE']=1.3

   def create_bubble(self):
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
      cb=val.reshape(2*ds+1,2*db+1,2*db+1)
      return cb

   def eighth_bubble(self):
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

   def whosmin(self,mat,ub,lb,delta):
     bord=np.array(self.energy.shape())
     ene=self.energy.data
     ub=np.array(ub)
     lb=np.array(lb)
     delta=np.array(delta)
     eub=ub + delta
     elb=lb + delta
     mub=ub-lb
     mlb=np.array([0,0,0])
     umask=eub > bord
     lmask=elb < 0
     mub[umask]-=eub[umask]-bord[umask]
     mlb[lmask]-=elb[lmask]
     eub[umask]=bord[umask]
     elb[lmask]=0
     eview=ene[elb[0]:eub[0],elb[1]:eub[1],elb[2]:eub[2]]
     mview=mat[mlb[0]:mub[0],mlb[1]:mub[1],mlb[2]:mub[2]]
     cmat=mview < eview
     eview[cmat]=mview[cmat]
     #self.energy.data[elb[0]:eub[0],elb[1]:eub[1],elb[2]:eub[2]]=eview
     
   def update_energies(self,lb,ub):
      mcb=self.data.get_slice(lb,ub)
      vv=self.eival
      ff=self.eifeat.T
      # this is one because we do not want to repeat the position (0,0,0)
      for i in range(1,vv.size):
         #print("."),
         sys.stdout.flush()
         mat=mcb/vv[i]
         delta=ff[i]
         # update in the eight directions
         # check limits on the current direction
         delta=np.array([1,1,1])*delta[0]
         self.whosmin(mat,ub,lb,delta)
         delta=np.array([1,1,-1])*delta[0]
         self.whosmin(mat,ub,lb,delta)
         delta=np.array([1,-1,1])*delta[0]
         self.whosmin(mat,ub,lb,delta)
         delta=np.array([1,-1,-1])*delta[0]
         self.whosmin(mat,ub,lb,delta)
         delta=np.array([-1,1,1])*delta[0]
         self.whosmin(mat,ub,lb,delta)
         delta=np.array([-1,1,-1])*delta[0]
         self.whosmin(mat,ub,lb,delta)
         delta=np.array([-1,-1,1])*delta[0]
         self.whosmin(mat,ub,lb,delta)
         delta=np.array([-1,-1,-1])*delta[0]
         self.whosmin(mat,ub,lb,delta)
      #print("DONE")

   def fit(self,cube,verbose=False,use_meta=True):
      
      # Set the RMS, or automatically find an estimate for it
      if not self.par.has_key('RMS'):
         rms=cube.estimate_rms()
         self.par['RMS']=rms
      plt.ion()
      # TODO: set parameters according to meta
      velres=self.par['VELORES']
      fwhmbeam=self.par['FWHMBEAM']
      bsize=self.par['BSIZE']
      nrms=self.par['NRMS']
      # Copy the supplied cube into a work cube which will hold the
      # residuals remaining after subtraction of the fitted Gaussians. 
      self.data=cube.copy()
      datamax=self.data.max()
      datamax=datamax[0]
      self.syn=cube.empty_like()
      self.energy=cube.copy()
      mas=np.isnan(self.energy.data)
      self.energy.data[np.logical_not(mas)]=datamax
      # Sigma values
      self.sb=bsize*fwhmbeam*self.FWHM_TO_SIGMA
      self.ss=bsize*velres*self.FWHM_TO_SIGMA
      # Deltas
      self.db=int(np.sqrt(-2*self.sb*self.sb*np.log(nrms*rms/datamax)))
      self.ds=int(np.sqrt(-2*self.ss*self.ss*np.log(nrms*rms/datamax)))
      log.info("Datamax =="+str(datamax))
      log.info("RMS ="+str(rms))
      log.info("Computed Deltas ="+str((self.db,self.ds)))
      (self.eival,self.eifeat)=self.eighth_bubble()
      lb=(0,0,0)
      ub=self.data.shape()
      self.update_energies(lb,ub)
      cb=self.create_bubble()
      plt.subplot(1,3,1)
      plt.imshow(self.data.get_stacked())
      plt.subplot(1,3,2)
      plt.imshow(self.energy.get_stacked())
      plt.subplot(1,3,3)
      plt.imshow(self.syn.get_stacked())
      plt.pause(1)
      delta=np.array([self.ds,self.db,self.db])
      iterate = True
      niter=0
      while iterate:
         # Report the iteration number to the user if required.
         niter+=1
         if verbose:
            log.info("Iteration: "+str(niter))
         y,xmax=self.energy.max()
         xmax=np.array(xmax)
         log.info("Maximum energy E = "+str(y)+" at "+str(xmax))
         rem=y - nrms*rms
         if rem <= nrms*rms:
            iterate=False
            break
         log.info("Remove E - "+str(nrms)+"*rms = "+str(rem)+" > "+str(nrms*rms))
         ub=xmax + delta + 1
         lb=xmax - delta
         #print ub - lb
         #print cb.shape

         self.data.add_flux(-rem*cb,lb,ub)
         self.syn.add_flux(rem*cb,lb,ub)
         log.info("Updating Energies")
         self.update_energies(lb,ub)
         plt.clf()
         plt.subplot(1,3,1)
         plt.imshow(self.data.get_stacked())
         plt.subplot(1,3,2)
         plt.imshow(self.energy.get_stacked())
         plt.subplot(1,3,3)
         plt.imshow(self.syn.get_stacked())
         plt.pause(0.00001)
 
      if verbose:
         log.info("PerfBubbleClump finished normally")



