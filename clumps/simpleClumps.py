import numpy as np
from collections import deque
from scipy.optimize import fmin_bfgs,check_grad,approx_fprime
#from scipy.optimize.linesearch import (line_search_BFGS, line_search_wolfe1, line_search_wolfe2, line_search_wolfe2 as line_search)
#from optimize import fmin_bfgs

import copy
import matplotlib.pyplot as plt
import sys
from astropy import log

K=4*np.log(2.0)


class SimpleClumps:

   def __init__(self):
      # Set the very 
      self.defaultParams()
   
   def defaultParams(self):
      self.par=dict()
      # Spectral Resoluion in pixels (smoothing function)
      self.par['VELORES']=2.0
      # Beam resoluion in pixels (smoothing function)
      self.par['FWHMBEAM']=2.0
      # The maximum allowed number of failed fits between succesful fits.
      self.par['MAXSKIP']=10
      # Maximum Clumps
      self.par['MAXCLUMPS']=sys.maxint
      # The iterative process ends when "npad" consecutive clumps all had peak
      # values below "peak_thresh" or all had areas below "area_thresh".
      self.par['NPAD']=10
      # The lower threshold for clump peaks to a user-specified multiple of the RMS noise.
      self.par['THRESH']=2.0
      # The lower threshold for clump area to a user-specified number of pixels.
      self.par['MINPIX']=3
      # The lowest value (normalised to the RMS noise level) at which
      # model Gaussians should be evaluated. 
      self.par['MODELMIN']=0.5
      # The max allowed fraction of bad pixels in a clump. 
      self.par['MAXBAD']=0.05
      # No.of standard deviations at which to reject peaks
      self.par['NSIGMA']=3.0
      # But reject peaks only if at least NPEAKS were found 
      self.par['NPEAKS']=9
      # Parameters which control the modification of the weights done by
      # the chi2 (this modification is meant to give low weights to pixels
      # which do not influence the Gaussian model
      self.par['NWF']= 10
      self.par['MINWF']=0.8
      self.par['MAXWF']=1.1
      # Maximum number of function evaluations to be used when fitting an
      # individual clump.
      self.par['MAXNF']=100
      # Chi-square stiffness parameter "Sa" which encourages the peak
      # amplitude of the fitted gaussian close to the maximum value in the
      # observed data.
      self.par['SA']=1.0
      # Chi-square stiffness parameter "Sb" which encourages the
      # background value to stay close to its initial value. This is an extra
      # stiffness added by DSB which is not in the Stutzki & Gusten paper. It
      # is used because the background value is usually determined by data
      # points which have very low weight and is thus poorly constrained. It
      # would thus be possibly to get completely erroneous background values
      # without this extra stiffness.
      self.par['SB']=0.1
      # Chi-square stiffness parameter "S0" which encourages the peak
      # amplitude of the fitted gaussian to be below the maximum value in the
      # observed data.
      self.par['S0']=1.0
      # Chi-square stiffness parameter "Sc" which encourages the peak
      # position of the fitted gaussian to be close to the peak position in the
      # observed data. 
      self.par['SC']=1.0
      # The ratio of the weighting function FWHM to the observed FWHM. */
      self.par['WWIDTH']=50.0
      # The value for which the weight is considered already zero 
      self.par['WMIN']=0.05


   def optimize(self):
      # Unpack used parameters
      wwidth=self.par['WWIDTH']
      velres=self.par['VELORES']
      beamfwhm=self.par['FWHMBEAM']

      fobs=np.zeros(3)
      fobs[0] = beamfwhm
      fobs[1] = beamfwhm
      fobs[2] = velres 

      beta=wwidth*fobs
      (ld,lu)=(np.rint(self.cval-beta),np.rint(self.cval+beta))
      lb=np.array([ld,lu]).min(axis=0)
      ub=np.array([ld,lu]).max(axis=0)
      lb=lb[::-1]
      ub=ub[::-1]
      print(lb,ub)
      self.val=self.data.get_slice(lb,ub).ravel()
      self.feat=self.data.get_index_features(lb,ub)
      
      xw_off=(self.feat[0] - self.cval[0])
      yw_off=(self.feat[1] - self.cval[1])
      vw_off=(self.feat[2] - self.cval[2])
      X=np.matrix([xw_off*self.val,yw_off*self.val,vw_off*self.val])
      print X
      S=X*X.T
      L=np.linalg.pinv(S)
      print(L,S)
      X=np.matrix([xw_off,yw_off,vw_off])
      print X.shape, L.shape
      B=(L*X).A*X.A
      B=B.sum(axis=0)
      print("B",B)
      E=np.exp(-K*B)
      m=(self.val/E).min()
      print m*E
      return m*E,lb,ub

           
   def fit(self,cube,verbose=False,use_meta=True):

      self.data=cube.copy()
      self.syn=cube.empty_like()
      niter=0
      while niter < 10:
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


