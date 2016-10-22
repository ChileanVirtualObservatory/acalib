import numpy as np
from collections import deque

from astropy.table import Table
from scipy.optimize import fmin_bfgs,check_grad,approx_fprime
#from scipy.optimize.linesearch import (line_search_BFGS, line_search_wolfe1, line_search_wolfe2, line_search_wolfe2 as line_search)
#from optimize import fmin_bfgs

import copy
import matplotlib.pyplot as plt
import sys
from astropy import log

K=4*np.log(2.0)


def jac_chi2(par,gc):
   # If the background is fixed, include zero background value
   val=np.nan*np.ones_like(par)
   if  par[ 3 ] > 0.0 and par[ 5 ] > 0.0 and par[ 8 ] > 0.0: 
      if gc.fixback:
         par=np.insert(par,1,0.0)
      # update computations if necesary
      gc.update_comp(par)
      val=gc.get_jaco(par)
   return val

def chi2(par,gc):
   val=1e1000
   # If the background is fixed, include zero background value
   if  par[ 3 ] > 0.0 and par[ 5 ] > 0.0 and par[ 8 ] > 0.0: 
      if gc.fixback:
         par=np.insert(par,1,0.0)
      # update computations if necesary
      gc.update_comp(par)
      val=gc.get_chi2(par)
   return val

class GaussClumps:

   def __init__(self):
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
      # self.par['MAXBAD']=0.05
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
      self.par['WWIDTH']=2.0
      # The value for which the weight is considered already zero 
      self.par['WMIN']=0.05

   def get_jaco(self,par):
     sa=self.par['SA']
     sb=self.par['SB']
     sc=self.par['SC']
     jaco=np.zeros(11)
     mod=np.zeros(11)
     t=self.peakfactor*self.expv
     jaco[0]=-2*t.dot(self.wres)/self.wsum
     jaco[1]=-2*self.wres.sum()/self.wsum
     ddx=self.X/self.sx2
     ddy=self.Y/self.sy2
     ddv=self.vt_off/self.sv2
     mterm=self.peak*self.expv
     t = -K*(-2*(ddx*self.cosv - ddy*self.sinv) + 2*par[9]*ddv)
     t *= mterm
     jaco[2]=-2*t.dot(self.wres)/self.wsum
     t = -K*(-2*ddx*ddx*par[3]) + self.f3
     t *= mterm
     jaco[3]=-2*t.dot(self.wres)/self.wsum
     t = -K*(-2*(ddx*self.sinv + ddy*self.cosv) + 2*par[10]*ddv)
     t *= mterm
     jaco[4]=-2*t.dot(self.wres)/self.wsum
     t = -K*(-2*ddy*ddy*par[5]) + self.f5
     t *= mterm
     jaco[5]=-2*t.dot(self.wres)/self.wsum
     t = -K*(-2*(ddx*(self.x_off*self.sinv - self.y_off*self.cosv) + ddy*(self.x_off*self.cosv + self.ym_off*self.sinv)))
     t *= mterm
     jaco[6]=-2*t.dot(self.wres)/self.wsum
     t = -K*(-2*ddv) 
     t *= mterm
     jaco[7]=-2*t.dot(self.wres)/self.wsum
     t = -K*(-2*ddv*ddv*par[8]) + self.f8
     t *= mterm
     jaco[8]=-2*t.dot(self.wres)/self.wsum
     t = -K*(-2*ddv*self.x_off)
     t *= mterm
     jaco[9]=-2*t.dot(self.wres)/self.wsum
     t = -K*(-2*ddv*self.y_off)
     t *= mterm
     jaco[10]=-2*t.dot(self.wres)/self.wsum
     # second pass
     jaco[0]+=2*sa*self.pdiff*self.peakfactor
     jaco[1]+=2*sa*self.pdiff + 2*sb*self.back_term
     jaco[2]+=2*4*sc*self.xm_off/self.bfsq
     jaco[3]+=2*sa*self.pdiff*self.f3*par[0]*self.peakfactor
     jaco[4]+=2*4*sc*self.ym_off/self.bfsq
     jaco[5]+=2*sa*self.pdiff*self.f5*par[0]*self.peakfactor
     jaco[7]+=2*4*sc*self.vm_off/self.velsq
     jaco[8]+=2*sa*self.pdiff*self.f8*par[0]*self.peakfactor
     if self.fixback:
        np.delete(jaco,[1])
     return jaco
   
   def get_chi2(self,par):
     sa=self.par['SA']
     sb=self.par['SB']
     sc=self.par['SC']
     chi2=self.wres.dot(self.res)
     chi2/=self.wsum
     off = (self.xm_off*self.xm_off + self.ym_off*self.ym_off )/self.bfsq
     off += self.vm_off*self.vm_off/self.velsq
     chi2 += sa*self.pdiff*self.pdiff + 4*sc*off +  sb*self.back_term*self.back_term
     return chi2
   
   def update_results(self,clump,lb,ub):
     self.update_comp(clump)
     ff=self.model - clump[1]
     #print lb,ub
     lb=self.data.fix_limits(lb)
     ub=self.data.fix_limits(ub)
     #print lb,ub
     #print ub[0]-lb[0],ub[1]-lb[1],ub[2]-lb[2]
     #print ff.shape
     ff=ff.reshape(ub[0]-lb[0],ub[1]-lb[1],ub[2]-lb[2])
     ff*=self.par['RMS']
     self.data.add_flux(-ff,lb,ub)
     ccode=self.caa.data.max() + 1
     caaff=self.caa.data[slice(lb[0],ub[0]),slice(lb[1],ub[1]),slice(lb[2],ub[2])]
     synff=self.syn.data[slice(lb[0],ub[0]),slice(lb[1],ub[1]),slice(lb[2],ub[2])]
     tmpff=ff.copy()
     tmpff[tmpff<self.par['RMS']]=0
     caaff[synff<tmpff]=ccode
     self.caa.replace_flux(caaff,lb,ub)
     self.syn.add_flux(ff,lb,ub)
     #TODO: improve area computation, right now using weigth ub, lb
     area=(ub-lb).sum()
     csum=ff.sum()
     return (csum,area)
   
   def update_comp(self,par):
     if np.array_equal(par,self.old_par):
         return
     self.old_par=par
     self.back_term=par[1] - self.guess[1]
     # Unpack parameters
     nwf=self.par['NWF']
     minwf=self.par['MINWF']
     maxwf=self.par['MAXWF']
     s0=self.par['S0']

     # Get the factor by which to correct the peak amplitude of the model to
     # take account of the smoothing by the instrumental beam.
     t = par[3]*par[3]
     sx2 = self.bfsq + t
     f3 = self.bfsq/(par[3]*sx2)
     peakfactor = t/sx2
     t = par[5]*par[5]
     sy2 = self.bfsq + t
     f5 = self.bfsq/(par[5]*sy2)
     peakfactor *= t/sy2
     t = par[8]*par[8]
     sv2 = self.velsq + t
     f8 = self.velsq/(par[8]*sv2)
     peakfactor *= t/sv2 
    
     if peakfactor > 0.0:
        peakfactor = np.sqrt(peakfactor)
     else:
        peakfactor = 0.0
     self.peak=par[0]*peakfactor
     self.sx2=sx2
     self.sy2=sy2
     self.sv2=sv2
     self.f3 = f3
     self.f5 = f5
     self.f8 = f8
     self.peakfactor=peakfactor

     # The difference between the model peak value (after being reduced to
     # take account of instrumental smoothing) and the data peak value.
     self.pdiff = self.peak + par[1] - self.valmax

     # The offset from the model centre to the data peak 
     xm_off = par[2] - self.cval[0]
     ym_off = par[4] - self.cval[1]
     vm_off = par[7] - self.cval[2]

     # Get the Gaussian model. Store the residual between the Gaussian model and data
     self.cosv = np.cos(par[6])
     self.sinv = np.sin(par[6])
     x_off=self.feat[0] - par[2]
     y_off=self.feat[1] - par[4]
     v_off=self.feat[2] - par[7]
     X = x_off*self.cosv + y_off*self.sinv
     Y = -x_off*self.sinv + y_off*self.cosv
     em = ( X*X/sx2 ) + ( Y*Y/sy2 )
     self.vt_off=v_off - par[9]*x_off - par[10]*y_off
     em += self.vt_off*self.vt_off/sv2
     expv = np.exp( -K*em )
     self.expv=expv
     model=self.peak*expv+ par[1]
     res= self.val - model
     self.X=X
     self.Y=Y
     self.x_off=x_off
     self.y_off=y_off
     self.v_off=v_off
     # If the changing of the model parameters make little difference to the
     # residuals at a given place in the data, then those residuals should be
     # given less weight since they could dominate the chi-squared value. If
     # the residual at the current pixel has not change by much since the
     # previous call, reduce the weight associated with the pixel. However,
     # if the parameter has not change by much then you would not expect the
     # residuals to change by much. Therefore, do not reduce the weight by so
     # much if the model value at this pixel has not changed by much since the
     # last call. In order to avoid instability, we only do this modification
     # for a few iterations near the start, and then allow the fitting
     # process to complete with fixed weights.
     if (not self.fixback) and (self.nf > 2) and (self.nwm <= nwf):
        # Only modify the weights if the background has changed. Without this,
        # the outlying background regions would be given low weights if the
        # background has not changed, resulting in the background being poorly
        # determined.
        if self.bg != 0.0:
           dbg=(par[1] - self.bg)/self.bg > 0.001
        else:
           dbg=(par[1] != 0.0)
        if dbg:
           wf=(res-self.res)/res
           wf/=(model - self.model)/model
           wf=np.abs(wf)
           wf[wf<minwf]=minwf
           wf[wf>maxwf]=maxwf
           wf[np.isnan(wf)]=1.0
           self.we*=wf
           self.we[self.we > 1.0]=1.0
           self.nwm+=1
     self.model=model 
     self.res=res
     self.xm_off=xm_off
     self.ym_off=ym_off
     self.vm_off=vm_off

     # Determine a scale factor which encourages the fitted intensity to stay
     # below the observed intensity. This does the same job as the
     # "s0.exp( Yi_fit - Yi )" term in the chi-squared expression given in
     # the Stutski & Gusten paper. The form used here was inherited from the
     # implementation of GaussClumps (obtained from
     # ftp.astro.uni-bonn.de/pub/heith/gaussclumps on 27/9/05) upon which this
     # implementation was based.
     rr = (s0+1)*np.ones_like(res)
     rr[res > 0.0]= 1.0

     # Compute the sum of chi-squared. We save the scaled residuals
     # in a work array (pr) so that we do not need to calculate them again if
     # this function is called subsequently to find the gradient for the same
     # set of parameer values.
     self.wsum = self.we.sum()
     self.wres = self.we*res*rr
     # Remember the background value for next time.
     self.bg = par[1]
     # Update nf
     self.nf+=1
   
   def optimize(self):
      # Unpack used parameters
      maxnf=self.par['MAXNF']
      wwidth=self.par['WWIDTH']
      wmin=self.par['WMIN']
      rms=self.par['RMS']
      velres=self.par['VELORES']
      beamfwhm=self.par['FWHMBEAM']
      self.bfsq=beamfwhm*beamfwhm
      self.velsq=velres*velres
      # Gaussian Window

      # The factor which scales the FWHM on each axis to the half-width of the
      # section of the data array to be be fitted. 
      beta=0.5*wwidth*np.sqrt(-np.log( wmin )/ np.log( 2.0 ) )
      (ld,lu)=(np.rint(self.cval-beta*self.fobs),np.rint(self.cval+beta*self.fobs))
      lb=np.array([ld,lu]).min(axis=0)
      ub=np.array([ld,lu]).max(axis=0)
      lb=lb[::-1]
      ub=ub[::-1]

      # Store the data normalised to the
      # RMS noise level. Also calculate and store the Gaussian weight for the
      # pixel.
      self.val=self.data.cut(lb,ub).copy().ravel()
      self.feat=self.data.index_features(lb,ub)
      
      xw_off=(self.feat[0] - self.cval[0])/(self.fobs[0]*wwidth)
      yw_off=(self.feat[1] - self.cval[1])/(self.fobs[1]*wwidth)
      vw_off=(self.feat[2] - self.cval[2])/(self.fobs[2]*wwidth)
      self.we=np.exp(-K*(xw_off*xw_off + yw_off*xw_off + vw_off*xw_off))
      self.we[self.we < wmin]=0.0
      
      # Normalise all other data values in the guess structure and in the 
      # array to the RMS noise level.
      self.val=self.val/rms
      self.valmax /= rms
      guess=self.guess
      guess[1] /= rms
      guess[0] /= rms

      # Number of invocations of the function
      self.nf=0
      self.nwm=0
      
      # Get the factor by which to correct the peak amplitude of the model to
      # take account of the smoothing by the instrumental beam.
      t = guess[3]*guess[3]
      dx0_sq = self.bfsq + t
      peakfactor = t/dx0_sq
      t = guess[5]*guess[5]
      dx1_sq = self.bfsq + t
      peakfactor *= t/dx1_sq
      t = guess[8]*guess[8]
      dv_sq = self.velsq + t
      peakfactor *= t/dv_sq

      # Do the correction.
      if  peakfactor > 0.0:
         guess[0] /= np.sqrt(peakfactor)
      
      if self.fixback:
         np.delete(guess,[1])

      marg=(self,)
      retval=fmin_bfgs(chi2, guess,args=marg,disp=False,full_output=True,fprime=jac_chi2)
      # Unpack results
      (xopt,fopt,gopt,Bopt,func_calls,grad_calls,warnflag)=retval
      if warnflag!=0 and self.fixback:
         self.fixback=False
         (xopt,fopt,gopt,Bopt,func_calls,grad_calls,warnflag)=fmin_bfgs(chi2, guess,args=marg,fprime=jac_chi2,maxiter=maxnf,disp=True)
      if self.fixback:
         np.insert(xopt,1,self.bg)
      if (xopt == self.guess).all():
         xopt=None
      # TODO: come back to normality!
      return xopt,lb,ub

   # TODO: Document this stuff (using cupid code...)
   def profWidth(self,dim):
      rms=self.par['RMS']
      if dim==0:
         vn=[0,0,-1]
         vp=[0,0,1]
         fwhm=self.par['FWHMBEAM']
      elif dim==1:
         vn=[0,-1,0]
         vp=[0,1,0]
         fwhm=self.par['FWHMBEAM']
      else:
         vn=[-1,0,0]
         vp=[1,0,0]
         fwhm=self.par['VELORES']

      # left search for significant minima
      left=np.array(self.imax)
      prev=np.nan
      vlow=self.data.data[self.imax]
      plow=self.imax
      csum=0.0
      nsum=0
      while True:
         left+=vn
         if (left < np.zeros(3)).any():
            left-=vn
            break
         val=self.data.data[tuple(left)]
         if np.ma.is_masked(val) or val==np.nan:
            prev=np.nan
            continue
         if val < vlow and prev != np.nan and prev - val < 1.5*rms:
            vlow=val
            plow=left.copy()
            #print "low cand",plow
            csum=0.0
            nsum=0
         else:
           csum+=val
           nsum+=1
           if csum/nsum - vlow >= 3*rms/np.sqrt(nsum) and nsum >= fwhm:
              break
         prev=val
      vlow+=rms
      #print "low",vlow,plow
      # Do the same working upwards from the peak to upper axis values.
      prev=np.nan
      vup=self.data.data[self.imax]
      pup=self.imax
      csum=0.0
      nsum=0
      right=np.array(self.imax)
      while True:
         right+=vp
         if (right >= np.array(self.data.data.shape)).any():
            right-=vp
            break
         val=self.data.data[tuple(right)]
         if np.ma.is_masked(val) or val==np.nan:
            prev=np.nan
            continue
         if val < vup and prev != np.nan and prev - val < 1.5*rms:
            vup=val
            pup=right.copy()
            csum=0.0
            nsum=0
         else:
           csum+=val
           nsum+=1
           if csum/nsum - vup >= 3*rms/np.sqrt(nsum) and nsum >= fwhm:
              break
         prev=val
      vup+=rms
      #print "up",vup,pup
      try:
         off=np.min(vlow,vup) + rms
      except ValueError:
        print(vlow,vup,rms)
        print(self.imax)
        self.data.data[self.imax]
        sys.exit()
      if vlow < vup:
         hgt=self.valmax - vlow
         cand=self.data.data[plow[0]:self.imax[0]+1,plow[1]:self.imax[1]+1,plow[2]:self.imax[2]+1]
         try:
           cand=cand[0][0]
         except IndexError:
           print(cand)
           print(plow, self.imax)
           sys.exit()
         cand-=vlow
         cand=cand[::-1]
         default=(self.imax-plow).sum()/2.0 
      else:
         hgt=self.valmax - vup
         cand=self.data.data[self.imax[0]:pup[0]+1,self.imax[1]:pup[1]+1,self.imax[2]:pup[2]+1]
         cand=cand[0][0]
         np.delete(cand,0)
         cand-=vup
         default=(np.array(self.imax)-np.array(pup)).sum()/2.0 
      cand=cand/hgt
      idx=np.arange(1,cand.size+1)
      cand=np.ma.fix_invalid(cand,fill_value=0.0)
      idx=idx[cand>0.25]
      cand=cand[cand>0.25]
      idx=idx[cand<0.75]
      cand=cand[cand<0.75]
      if cand.size!=0:
         default=1.665*(idx/np.log(cand)).sum()/cand.size
      return (default,off)
      
   def setInit(self):
      # Unpack used parameters
      beamfwhm=self.par['FWHMBEAM']
      velres=self.par['VELORES']
      rms=self.par['RMS']

      guess=np.zeros(11)      
      # Get a guess at the observed clump fwhm by forming a radial profile and
      # finding the distance to the first significant minimum. This also increments
      # "off" by the minimum (i.e. base line) data value in the profile. Do
      # this for both spatial axes, and then take the mean (i.e. we assume the
      # clump is circular as an initial guess)
      self.fobs=np.zeros(3)
      off=np.zeros(3)
      self.fobs[0],off[0] = self.profWidth(0) 
      self.fobs[1],off[1] = self.profWidth(1) 
      self.fobs[2],off[2] = self.profWidth(2)
      # TODO: Small Hack
      if self.fobs[2]<velres+0.1:
          self.fobs[2]=velres+0.1
      fbeam=0.5*(self.fobs[0]  + self.fobs[1])/beamfwhm
      if fbeam < 1.0: 
         fbeam=1.2
      self.fobs[0] = fbeam*beamfwhm
      self.fobs[1] = fbeam*beamfwhm
           
      # Store the Guessed model
      self.cval=np.array([self.imax[2],self.imax[1],self.imax[0]])
      guess[2]=self.cval[0]
      guess[4]=self.cval[1]
      guess[7]=self.cval[2]
      # Find the initial guess at the intrinsic FWHM (i.e. the FWHM of the
      # clump before being blurred by the instrument beam). Do the same for 
      # the second axis. Assume zero rotation of the elliptical clump shape.
      guess[3]=np.sqrt(fbeam*fbeam- 1.0 )*beamfwhm
      guess[5]=guess[3]
      guess[6]=0.0
      # Now do the same for the third (velocity) axis if necessary. Assume
      # zero velocity gradient
      fvel=self.fobs[2]/velres
      guess[8]=np.sqrt(fvel*fvel- 1.0 )*velres
      if np.isnan(guess[8]):
         print(fvel,self.fobs[2],fvel*fvel)
      guess[9]=0.0
      guess[10]=0.0
      
      # Store the mean of the background estimates, and the peak value. Noise
      # will result in the peak data value being larger than the peak clump value
      # by about the RMS noise. Therefore, reduce the peak value by the RMS.
      guess[1] = off.sum()/3
      guess[0] = self.valmax - guess[1] - rms

      # Negative background levels are unphysical (since it is assumed that
      # any background has already been removed from the data before running
      # this algorithm (TODO: CHECK)). However, an apparent negative background can be formed by
      # a previous ill-position fit resulting in negative residiauls. Therefore
      # we have to guard against negative backgrounds. If the initial background
      # estimate is significantly less than zero, then set it to zero, and
      # indicate that the background value should be fixed (i.e. not included
      # as a free parameter in the fitting process). Here, "significant" means
      # more than 5% of the total peak height. */
      self.fixback=False
      if guess[1] < -np.abs(guess[ 0 ]*0.05):
         guess[0] += guess[1]
         guess[1] = 0.0
         self.fixback = True
      self.guess=guess
      self.old_par=None

   def fit(self,cube,verbose=False,use_meta=True):

      FWHM_TO_SIGMA = 1. / (8 * np.log(2))**0.5
      # Set the RMS, or automatically find an estimate for it
      if not self.par.has_key('RMS'):
         rms=cube.rms()
         self.par['RMS']=rms
      
      # TODO: set parameters according to meta
 
      # Unpack used parameters
      npeaks=self.par['NPEAKS']
      mlim=self.par['MODELMIN']
      peak_thresh=self.par['THRESH']
      area_thresh=self.par['MINPIX']
      maxclump=self.par['MAXCLUMPS']
      npad=self.par['NPAD']
      maxskip=self.par['MAXSKIP']
      nsig=self.par['NSIGMA']

      # Copy the supplied cube into a work cube which will hold the
      # residuals remaining after subtraction of the fitted Gaussians. 
      self.data=cube.copy()
      self.syn=cube.empty_like()
      self.caa=cube.empty_like()
      # Initialise the number of clumps found so far.
      iclump = 0

      # Indicate that no peaks have been found below the lower threshold for clump
      # peak values, or below the lower area threshold. 
      peaks_below = 0
      area_below = 0
 
      # Initialise the variables used to keep track of the mean and standard
      # deviation of the most recent "npeak" fitted peak values. 
      mean_peak = 0.0
      sigma_peak = 0.0 
      # The value most recently added to "peaks"
      new_peak = 0.0
      # Sum of the values in "peaks"
      sum_peak = 0.0
      # Sum of the squares of the values in "peaks" 
      sum_peak2 = 0.0

      # Number of pixels contributing to the clump
      area=0

      # Iterations performed so far 
      niter = 0
      iterate = True
      # No. of failed fits since last good fit 
      nskip = 0
      # Sum of the values in all the used clumps so far 
      sumclumps = 0.0
      # Sum of the supplied data values 
      sumdata = self.data.flux()
      
      # peaks contains the last npeaks... 
      peaks=np.zeros(npeaks)
      clist=Table(names=("Intensity", "Offset", "RA mu", "RA std","DEC mu", "DEC std","Angle","FREQ mu","FREQ std","RA vel grad","DEC vel grad"),dtype=('f8','f8','f8','f8','f8','f8','f8','f8','f8','f8','f8'))
      clist.meta=cube.meta
      # Loop round fitting a gaussian to the largest remaining peak in the
      # residuals array. */
      while iterate:
         # Report the iteration number to the user if required.
         niter+=1
         if verbose:
            log.info("Iteration: "+str(niter))
         # Find the cube index of the element with the largest value in the residuals cube.
         # imax: Index of element with largest residual
         (self.valmax,self.imax) = self.data.max()
         # Finish iterating if all the residuals are bad, or if too many iterations
         # have been performed since the last succesfully fitted clump. 
         if np.isnan(self.imax).any():
            iterate = False
            niter-=1
            if verbose:
               log.info("There are no good pixels left to be fitted.")
               continue
         elif  nskip > maxskip:
            iterate = False
            niter-=1
            if verbose:
               log.info("The previous "+str(maxskip)+" fits were unusable.")
               continue
         # If not, make an initial guess at the Gaussian clump parameters centred on the current peak.
         self.setInit()

         
         # Find the best fitting parameters, starting from the above initial guess.
         (clump,lb,ub)=self.optimize()
         # If no fit could be performed, then found = False
         if clump!=None:
            # Skip this fit if we have an estimate of the standard deviation of the
            # "npeaks" most recent clump peak values, and the peak value of the clump
            # just fitted is a long way (more than NSIGMA standard deviations) from the
            # peak value of the previously fitted clump. Also skip it if the peak
            # value is less than the "mlim" value.
            if (peaks.size == 0 or iclump < npeaks or np.abs(clump[0] - new_peak) < nsig*sigma_peak ) and clump[0] > mlim:

               # Record the new peak value for use with the next peak, and update the
               # standard deviation of the "npeaks" most recent peaks. These values are
               # stored cyclically in the "peaks" array. */
               if peaks.size > 0:
                  np.roll(peaks,1)
                  new_peak = clump[0]
                  old_peak = peaks[0]
                  peaks[0]=new_peak
                  sum_peak += new_peak - old_peak
                  sum_peak2 += new_peak*new_peak - old_peak*old_peak
                  if sum_peak2 < 0.0:
                     sum_peak2 = 0.0
                  mean_peak = sum_peak/npeaks
                  sigma_peak = np.sqrt(sum_peak2/npeaks - mean_peak*mean_peak)
               
               # Increment the number of peaks found. 
               iclump+=1

               # Reset the number of failed fits since the last good fit. */
               nskip = 0

               # Remove the model fit (excluding the background) from the residuals.
               # This also creates data values asociated with the clumps
               # The standard deviation of the new residuals is returned. */
               if clump[0] >= peak_thresh:
                  #record clump
                  clist+=tuple(clump)
                  (csum,area)=self.update_results(clump,lb,ub)
                  sumclumps+=csum
                  # TODO: implement this!
                  # Display the clump parameters on the screen if required. */
                  #cupidGCListClump( iclump, ndim, x, chisq, slbnd, rms, status )

               # If this clump has a peak value which is below the threshold, increment
               # the count of consecutive clumps with peak value below the threshold.
               # Otherwise, reset this count to zero.
               if clump[0] < peak_thresh:
                  self.data.data.mask[self.imax]=True
                  peaks_below+=1
               else:
                  peaks_below=0

               # If this clump has an area which is below the threshold, increment
               # the count of consecutive clumps with area below the threshold.
               # Otherwise, reset this count to zero. 
               if area < area_thresh:
                  area_below+=1
               else:
                  area_below=0

               # If the maximum number of clumps have now been found, exit.*/
               if iclump == maxclump:
                  iterate = False
                  if verbose:
                     log.info("The specified maximum number of clumps ("+str(maxclump)+") have been found.")

               # If the integrated data sum in the fitted gaussians exceeds or equals
               # the integrated data sum in the input, exit. 
               elif sumclumps >= sumdata:
                  iterate = False
                  if verbose:
                     log.info("The total data sum of the fitted Gaussians ("+str(sumclumps)+") has reached the total data sum in the supplied data ("+str(sumdata)+").")
 
               # If the count of consecutive peaks below the threshold has reached
               # "Npad", terminate.
               elif peaks_below == npad:
                  iterate = False
                  if verbose:
                     log.info("The previous"+str(npad)+"clumps all had peak values below the threshold.")

               # If the count of consecutive clumps with area below the threshold has reached
               # "Npad", terminate.
               elif area_below == npad:
                  iterate = False
                  if verbose:
                     log.info("The previous "+str(npad)+" clumps all had areas below the threshold.")
               
            # If the peak value fitted is very different from the previous fitted peak
            # value, set the residuals array element bad in order to prevent the
            # algorithm from trying to fit a peak to the same pixel again. 
            else:
               self.data.data.mask[self.imax]=True
               new_peak = 0.5*(new_peak + clump[0])
               nskip+=1
               if verbose:
                  log.info("Clump rejected due to aberrant peak value. Ignoring Pixel...")

         # Tell the user if no clump could be fitted around the current peak
         # pixel value 
         else:
            nskip+=1
            if verbose:
               log.info("No clump fitted (optimization falied). Ignoring Pixel...")
            # Set the specified element of the residuals array bad if no fit was
            # performed. This prevents the any subsequent attempt to fit a Gaussian
            # to the same peak value.
            self.data.data.mask[self.imax]=True
      if verbose:
         log.info("GaussClump finished normally")
         # TODO: Usable Clumps
         ## Tell the problems with the clumps. */
         #if nclump == 0:
         #   print "No usable clumps found."
         #if iclump - nclump >= 1:
         #   print iclump - nclump,"clump(s) rejected because they touch an edge of the data array."
         ## Tell the user how many iterations have been performed (i.e. how many
         ## attempts there have been to fit a Gaussian peak
         #if niter == 1:
         #   print "No fit attempted."
         #else:
         #   print "Fits attempted for ",iclump," candidate clumps (",niter-iclump," failed)."
      return clist
        


#   def fit_up_to_SNR(self,cube,snr,verbose=False,):
#
#      FWHM_TO_SIGMA = 1. / (8 * np.log(2))**0.5
#      # Set the RMS, or automatically find an estimate for it
#      if not self.par.has_key('RMS'):
#         rms=cube.estimate_rms()
#         self.par['RMS']=rms
#      
#      # TODO: set parameters according to meta
# 
#      # Unpack used parameters
#      npeaks=self.par['NPEAKS']
#      mlim=self.par['MODELMIN']
#      peak_thresh=self.par['THRESH']
#      area_thresh=self.par['MINPIX']
#      maxclump=self.par['MAXCLUMPS']
#      npad=self.par['NPAD']
#      maxskip=self.par['MAXSKIP']
#      nsig=self.par['NSIGMA']
#
#      # Copy the supplied cube into a work cube which will hold the
#      # residuals remaining after subtraction of the fitted Gaussians. 
#      self.data=cube.copy()
#      self.syn=cube.empty_like()
#      # Initialise the number of clumps found so far.
#      iclump = 0
#
#      # Indicate that no peaks have been found below the lower threshold for clump
#      # peak values, or below the lower area threshold. 
#      peaks_below = 0
#      area_below = 0
# 
#      # Initialise the variables used to keep track of the mean and standard
#      # deviation of the most recent "npeak" fitted peak values. 
#      mean_peak = 0.0
#      sigma_peak = 0.0 
#      # The value most recently added to "peaks"
#      new_peak = 0.0
#      # Sum of the values in "peaks"
#      sum_peak = 0.0
#      # Sum of the squares of the values in "peaks" 
#      sum_peak2 = 0.0
#
#      # Number of pixels contributing to the clump
#      area=0
#
#      # Iterations performed so far 
#      niter = 0
#      iterate = True
#      # No. of failed fits since last good fit 
#      nskip = 0
#      # Sum of the values in all the used clumps so far 
#      sumclumps = 0.0
#      # Sum of the supplied data values 
#      sumdata = self.data.flux()
#      
#      # peaks contains the last npeaks... 
#      peaks=np.zeros(npeaks)
#      clist=[]
#      # Loop round fitting a gaussian to the largest remaining peak in the
#      # residuals array. */
#      while iterate:
#         # Report the iteration number to the user if required.
#         niter+=1
#         if verbose:
#            log.info("Iteration: "+str(niter))
#         # Find the cube index of the element with the largest value in the residuals cube.
#         # imax: Index of element with largest residual
#         (self.valmax,self.imax) = self.data.max()
#         if verbose:
#              print "GAP ",self.valmax - (rms + snr*rms)
#         if self.valmax < rms   + snr*rms:
#            iterate = False
#            niter-=1
#            if verbose:
#               log.info("Maximum point is at SNR limit.")
#         # Finish iterating if all the residuals are bad, 
#         if np.isnan(self.imax).any():
#            iterate = False
#            niter-=1
#            if verbose:
#               log.info("There are no good pixels left to be fitted.")
#               continue
#         # If not, make an initial guess at the Gaussian clump parameters centred on the current peak.
#         self.setInit()
#         
#         # Find the best fitting parameters, starting from the above initial guess.
#         (clump,lb,ub)=self.optimize()
#         # If no fit could be performed, then found = False
#         if clump!=None:
#            # Skip this fit if we have an estimate of the standard deviation of the
#            # "npeaks" most recent clump peak values, and the peak value of the clump
#            # just fitted is a long way (more than NSIGMA standard deviations) from the
#            # peak value of the previously fitted clump. Also skip it if the peak
#            # value is less than the "mlim" value.
#            if (peaks.size == 0 or iclump < npeaks or np.abs(clump[0] - new_peak) < nsig*sigma_peak ) and clump[0] > mlim:
#
#               # Record the new peak value for use with the next peak, and update the
#               # standard deviation of the "npeaks" most recent peaks. These values are
#               # stored cyclically in the "peaks" array. */
#               if peaks.size > 0:
#                  np.roll(peaks,1)
#                  new_peak = clump[0]
#                  old_peak = peaks[0]
#                  peaks[0]=new_peak
#                  sum_peak += new_peak - old_peak
#                  sum_peak2 += new_peak*new_peak - old_peak*old_peak
#                  if sum_peak2 < 0.0:
#                     sum_peak2 = 0.0
#                  mean_peak = sum_peak/npeaks
#                  sigma_peak = np.sqrt(sum_peak2/npeaks - mean_peak*mean_peak)
#               
#               # Increment the number of peaks found. 
#               iclump+=1
#
#               # Reset the number of failed fits since the last good fit. */
#               nskip = 0
#
#               # Remove the model fit (excluding the background) from the residuals.
#               # This also creates data values asociated with the clumps
#               # The standard deviation of the new residuals is returned. */
#               if clump[0] >= peak_thresh:
#                  #record clump
#                  clist.append(clump)
#                  (csum,area)=self.updateResults(clump,lb,ub)
#                  sumclumps+=csum
#                  # TODO: implement this!
#                  # Display the clump parameters on the screen if required. */
#                  #cupidGCListClump( iclump, ndim, x, chisq, slbnd, rms, status )
#
#               # If this clump has a peak value which is below the threshold, increment
#               # the count of consecutive clumps with peak value below the threshold.
#               # Otherwise, reset this count to zero.
#               if clump[0] < peak_thresh:
#                  self.data.data.mask[self.imax]=True
#                  if verbose:
#                     log.info("Clump rejected because it is below threshold.")
#                  peaks_below+=1
#               else:
#                  peaks_below=0
#
#               # If this clump has an area which is below the threshold, increment
#               # the count of consecutive clumps with area below the threshold.
#               # Otherwise, reset this count to zero. 
#               if area < area_thresh:
#                  log.info("Clump rejected because it is below threshold.")
#                  area_below+=1
#               else:
#                  area_below=0
#
#               # If the maximum number of clumps have now been found, exit.*/
#               if iclump == maxclump:
#                  iterate = False
#                  if verbose:
#                     log.info("The specified maximum number of clumps ("+str(maxclump)+") have been found.")
#
#               # If the integrated data sum in the fitted gaussians exceeds or equals
#               # the integrated data sum in the input, exit. 
#               elif sumclumps >= sumdata:
#                  iterate = False
#                  if verbose:
#                     log.info("The total data sum of the fitted Gaussians ("+str(sumclumps)+") has reached the total data sum in the supplied data ("+str(sumdata)+").")
# 
#            # If the peak value fitted is very different from the previous fitted peak
#            # value, set the residuals array element bad in order to prevent the
#            # algorithm from trying to fit a peak to the same pixel again. 
#            else:
#               self.data.data.mask[self.imax]=True
#               new_peak = 0.5*(new_peak + clump[0])
#               nskip+=1
#               if verbose:
#                  log.info("Clump rejected due to aberrant peak value. Ignoring Pixel...")
#
#         # Tell the user if no clump could be fitted around the current peak
#         # pixel value 
#         else:
#            nskip+=1
#            if verbose:
#               log.info("No clump fitted (optimization falied). Ignoring Pixel...")
#            # Set the specified element of the residuals array bad if no fit was
#            # performed. This prevents the any subsequent attempt to fit a Gaussian
#            # to the same peak value.
#            self.data.data.mask[self.imax]=True
#      if verbose:
#         log.info("GaussClump finished normally")
#         # TODO: Usable Clumps
#         ## Tell the problems with the clumps. */
#         #if nclump == 0:
#         #   print "No usable clumps found."
#         #if iclump - nclump >= 1:
#         #   print iclump - nclump,"clump(s) rejected because they touch an edge of the data array."
#         ## Tell the user how many iterations have been performed (i.e. how many
#         ## attempts there have been to fit a Gaussian peak
#         #if niter == 1:
#         #   print "No fit attempted."
#         #else:
#         #   print "Fits attempted for ",iclump," candidate clumps (",niter-iclump," failed)."
#      return clist

