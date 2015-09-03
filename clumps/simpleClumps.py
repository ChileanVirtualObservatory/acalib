import numpy as np
from collections import deque
from scipy.optimize import fmin_bfgs,check_grad,approx_fprime
#from scipy.optimize.linesearch import (line_search_BFGS, line_search_wolfe1, line_search_wolfe2, line_search_wolfe2 as line_search)
#from optimize import fmin_bfgs
import time
import copy
import matplotlib.pyplot as plt
import sys
from astropy import log
from mpl_toolkits.mplot3d import axes3d

K=4*np.log(2.0)


class SimpleClumps:

   def __init__(self):
      # Set the very 
      self.defaultParams()
   
   def defaultParams(self):
      plt.ion()
      self.par=dict()
      # Spectral Resoluion in pixels (smoothing function)
      self.par['VELORES']=2.0
      # Beam resoluion in pixels (smoothing function)
      self.par['FWHMBEAM']=4.0
      # The ratio of the weighting function FWHM to the observed FWHM. */
      self.par['WWIDTH']=5.0


   def optimize(self):
      # Unpack used parameters
      wwidth=self.par['WWIDTH']
      velres=self.par['VELORES']
      beamfwhm=self.par['FWHMBEAM']
      rms=self.par['RMS']

      fobs=np.zeros(3)
      fobs[0] = beamfwhm
      fobs[1] = beamfwhm
      fobs[2] = velres 

      beta=wwidth*fobs
      (ld,lu)=(np.rint(self.cval-beta),np.rint(self.cval+beta))
      lb=np.array([ld,lu]).min(axis=0)
      ub=np.array([ld,lu]).max(axis=0)
      lb=lb[::-1]
      ub=ub[::-1] +  1
      #print(lb,ub)
      #print self.imax
      cc=self.data.get_slice(lb,ub)
      self.val=cc.copy().ravel()
      self.feat=self.data.get_index_features(lb,ub)
      xw_off=(self.feat[0] - self.cval[0])
      yw_off=(self.feat[1] - self.cval[1])
      vw_off=(self.feat[2] - self.cval[2])
      X=np.matrix([xw_off,yw_off,vw_off])
      dist=(X.A*X.A).sum(axis=0)
      order=np.argsort(dist)
      
      for elm in order:
         if dist[elm]==0.0:
            continue
         vect=[-vw_off[elm],-yw_off[elm],-xw_off[elm]]
         fact=np.abs(vect)
         fact/=fact.sum()
         vect=np.sign(vect)
         vvv=np.unravel_index(elm,cc.shape)
         val=fact[0]*cc[vvv[0]+vect[0],vvv[1],vvv[2]] + fact[1]*cc[vvv[0],vvv[1]+vect[1],vvv[2]] + fact[2]*cc[vvv[0],vvv[1],vvv[2]+vect[2]]
         cc[vvv[0],vvv[1],vvv[2]]=min(cc[vvv[0],vvv[1],vvv[2]],val)
   
      ccc=cc.ravel()
      ccc[ccc<0.3*self.valmax]=0
      #plt.plot(dist)
      #plt.show()
      #plt.pause(10)
      #print xw_off.shape,yw_off.shape,vw_off.shape
      #print self.val.shape
      #plt.plot(vw_off)
      #plt.show()
      #plt.pause(10)
      #vv=vv*np.exp(-0.5*(xw_off*xw_off/xw_off.max() + yw_off*yw_off/yw_off.max() + vw_off*vw_off/vw_off.max()))
      #print vv.max()
      sqvv=np.sqrt(ccc)
      XP=np.matrix([sqvv*xw_off,sqvv*yw_off,sqvv*vw_off])
      S=XP*XP.T/np.square(sqvv).sum()
      #print "S",S
      L=np.linalg.pinv(S)
      #print "L",L
      #print X.shape, L.shape
      B=(L*X).A*X.A
      B=B.sum(axis=0)
      E=np.exp(-0.5*B)
      #plt.figure()
      #plt.imshow(PP.sum(axis=0))
      control=(E > 5*rms)
      m=np.nanmin((self.val[control])/E[control])
      print "rem", m 
      print "max", self.valmax
      res=m*E
      print "nonzero",100.0*(res>0).sum()/res.size,"%"
      res[np.logical_not(control)]=0.0
      #PP=res.reshape(ub[0]-lb[0],ub[1]-lb[1],ub[2]-lb[2])
      plt.clf()
      ax=plt.subplot(1, 3, 1,projection='3d')
      #OO=self.val.reshape(ub[0]-lb[0],ub[1]-lb[1],ub[2]-lb[2])
      #plt.imshow(OO.sum(axis=0))
      #plt.imshow(self.data.get_slice(lb,ub).sum(axis=0))
      control=(self.val>0.5*self.valmax)
      ax.set_xlim([self.feat[0].min(),self.feat[0].max()])
      ax.set_ylim([self.feat[1].min(),self.feat[1].max()])
      ax.set_zlim([self.feat[2].min(),self.feat[2].max()])
      ax.scatter(self.feat[0][control], self.feat[1][control], self.feat[2][control], c=self.val[control],marker='+')
      cm=plt.get_cmap()
      ax=plt.subplot(1, 3, 2,projection='3d')
      #control=(res>0.5*res.max())
      control=(ccc>0.0)
      ax.set_xlim([self.feat[0].min(),self.feat[0].max()])
      ax.set_ylim([self.feat[1].min(),self.feat[1].max()])
      ax.set_zlim([self.feat[2].min(),self.feat[2].max()])
      #ax.scatter(self.feat[0][control], self.feat[1][control], self.feat[2][control], c=res[control],marker='+',cmap=cm)
      ax.scatter(self.feat[0][control], self.feat[1][control], self.feat[2][control], c=ccc[control],marker='+',cmap=cm)
      #plt.imshow(PP.sum(axis=0))
      ax=plt.subplot(1, 3, 3,projection='3d')
      control=(res>0.0)
      ax.set_xlim([self.feat[0].min(),self.feat[0].max()])
      ax.set_ylim([self.feat[1].min(),self.feat[1].max()])
      ax.set_zlim([self.feat[2].min(),self.feat[2].max()])
      ax.scatter(self.feat[0][control], self.feat[1][control], self.feat[2][control], c=res[control],marker='+',cmap=cm)
      plt.show()
      plt.pause(0.001)
      #self.val[self.val==0]=np.nan
      #NN=E[self.val<0.2*self.valmax]
      #if (NN==np.nan).all():
      #   print "WTF!"
      #   return 0.0*E,lb,ub
      #print NN
      #tabu=np.nanmax(NN)
      #E[E<tabu]=np.nan
      #m=np.nanmin(self.val/E)
      #if m==0:
      #   print E[np.logical_not(np.isnan(E))]
      #print "imax",self.imax
      #print "tabu",tabu
      #print "m",m 
      #E[np.isnan(E)]=0.0
      #print m
      #plt.figure()
      #plt.plot(self.val,'r')
      #D=E.copy()
      #D[np.isnan(D)]=0.0
      #plt.plot(self.val - m*D ,'b')
      #plt.show()
      return res,lb,ub

           
   def fit(self,cube,verbose=False,use_meta=True):

      self.data=cube.copy()
      self.syn=cube.empty_like()
      if not self.par.has_key('RMS'):
         rms=cube.estimate_rms()
         self.par['RMS']=rms
      #dat=self.data.data.copy()
      #dat=dat - dat.min()
      #dat[dat<10*rms]=0
      #dat[dat>=3*rms]=1
      #plt.imshow(dat.sum(axis=0))
      #plt.show()
      niter=0
      while niter < 50:
         # Report the iteration number to the user if required.
         niter+=1
         if verbose:
            log.info("Iteration: "+str(niter))
         # Find the cube index of the element with the largest value in the residuals cube.
         # imax: Index of element with largest residual
         (self.valmax,self.imax) = self.data.max()
         self.cval=np.array([self.imax[2],self.imax[1],self.imax[0]])
         (ff,lb,ub)=self.optimize()
         ff=ff.reshape(ub[0]-lb[0],ub[1]-lb[1],ub[2]-lb[2])
         self.data.add_flux(-ff,lb,ub)
         self.syn.add_flux(ff,lb,ub)
      plt.ioff()


