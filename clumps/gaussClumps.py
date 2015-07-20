import numpy as np
from collections import deque
from scipy.optimize import fmin_bfgs,check_grad
import copy
import matplotlib.pyplot as plt
import sys
from cube import *

def chi2(par,gc):
   ret=None
   # If the background is fixed, include zero background value
   if gc.fixback:
      par=np.insert(par,1,0.0)
      back_term=0
   else
      back_term=par[1] - gc.guess[1]
      back_term *= gc.pars['SB']*back_term
   if not gc.reuse:
      if  par[ 3 ] <= 0.0: return ret
      if  par[ 5 ] <= 0.0: return ret
      if  par[ 8 ] <= 0.0: return ret
 
      # Get the factor by which to correct the peak amplitude of the model to
      # take account of the smoothing by the instrumental beam.
      t = par[ 3 ]*par[ 3 ]
      dx_sq = gc.bfsq + t
      peakfactor = t/dx_sq
      f3 = par[ 0 ]*gc.bfsq/( par[ 3 ]*dx_sq )
      t = par[ 5 ]*par[ 5 ]
      dx_sq = gc.bfsq + t
      peakfactor *= t/dx_sq;
      f5 = par[ 0 ]*gc.bfsq/( par[ 5 ]*dx_sq )
      t = par[ 8 ]*par[ 8 ]
      dx_sq = gc.velsq + t;
      peakfactor *= t/dx_sq;
      f8 = par[ 0 ]*gc.velsq/( par[ 8 ]*dx_sq );

      if peakfactor > 0.0:
         peakfactor = np.sqrt( peakfactor )
      else:
         peakfactor = 0.0;
      }

      f3 *= peakfactor;
      f5 *= peakfactor;
      f8 *= peakfactor;
 
      # The difference between the model peak value (after being reduced to
      # take account of instrumental smoothing) and the data peak value.
      pdiff = peakfactor*par[ 0 ] + par[ 1 ] - gc.ymax;
 
      # The offset from the model centre to the data peak 
      x0_off = par[ 2 ] - gc.cval[ 0 ];
      x1_off = par[ 4 ] - gc.cval[ 1 ];
      v_off = par[ 7 ] - gc.cval[ 2 ];

      # Initialise the total chi squared value */
      chisq = 0.0;

      # Get the Gaussian model. Store
      # the residual between the Gaussian model and data
      
292          m = cupidGCModel( ndim, x, par, -1, 1, ( iel == 0 ), status );
293          res = *py - m;
294 
295 /* If the changing of the model parameters make little difference to the
296    residuals at a given place in the data, then those residuals should be
297    given less weight since they could dominate the chi-squared value. If
298    the residual at the current pixel has not change by much since the
299    previous call, reduce the weight associated with the pixel. However,
300    if the parameter has not change by much then you would not expect the
301    residuals to change by much. Therefore, do not reduce the weight by so
302    much if the model value at this pixel has not changed by much since the
303    last call. In order to avoid instability, we only do this modification
304    for a few iterations near the start, and then allow the fitting
305    process to complete with fixed weights. */
306          if( !cupidGC.fixback && cupidGC.nf > 2 && nwm <= cupidGC.nwf ) {
307             if( res != 0.0 && m != 0.0 && m != *pm ) {
308 
309 /* Only modify the weights if the background has changed. Without this,
310    the outlying background regions would be given low weights if the
311    background has not changed, resulting in the background being poorly
312    determined. */
313                if( bg != 0.0 ) {
314                   dbg = ( fabs( ( par[ 1 ] - bg )/bg ) > 0.001 );
315                } else {
316                   dbg = ( par[ 1 ] != 0.0 );
317                }
318                if( dbg ) {
319                   wf = ( res - *pu )/ res;
320                   wf /= ( m - *pm )/ m;
321                   wf = fabs( wf );
322                   wf = ( wf < cupidGC.minwf ) ? cupidGC.minwf : ( wf > cupidGC.maxwf ) ? cupidGC.maxwf : wf ;
323                   *pw *= wf;
324                   if( *pw > 1.0 ) *pw = 1.0;
325                   wmod = 1;
326                }
327             }
328          }
329 
330 /* Save the residual and model value at this pixel */
331          *pu = res;
332          *pm = m;
 /* Determine a scale factor which encourages the fitted intensity to stay
335    below the observed intensity. This does the same job as the
336    "s0.exp( Yi_fit - Yi )" term in the chi-squared expression given in
337    the Stutski & Gusten paper. The form used here was inherited from the
338    implementation of GaussClumps (obtained from
339    ftp.astro.uni-bonn.de/pub/heith/gaussclumps on 27/9/05) upon which this
340    implementation was based. */
341          rr = ( res > 0.0 ) ? 1.0 : cupidGC.s0p1;
342 
343 /* Increment the running sum of chi-squared. We save the scaled residuals
344    in a work array (pr) so that we do not need to calculate them again if
345    this function is called subsequently to find the gradient for the same
346    set of parameer values. */
347          wsum += *pw;
348          *pr = *pw*res*rr;
349          chisq += *pr*res;
350          *prs = *pr*res;
351 
352 /* Move the pointers on to the next pixel in the section of the data
353    array being fitted. */
354          py++;
355          pw++;
356          pr++;
357          pu++;
358          pm++;
359          prs++;
360 
361 /* Get the grid coords (within the full size original data array) of the
362    next pixel in the section currently being fitted. This assumes fortran
363    ordering of the elements in the arrays.*/
364          iax = 0;
365          x[ iax ] += 1.0;
366          while( x[ iax ] > cupidGC.ubnd[ iax ] ) {
367             x[ iax ] = cupidGC.lbnd[ iax ];
368             if( ++iax == ndim ) break;
369             x[ iax ] += 1.0;
370          }
371       }
372 
373 /* Remember the background value for next time. */
374       bg = par[ 1 ];
375 
376 /* Increment the number of iteration sthat have made modifications to the
377    weights (if any such change has in fact been made). */
378       if( wmod ) nwm++;
379 
380 /* Divide by the sum of the weights . */
381       cupidGC.wsum = wsum;
382       chisq /= wsum;
384 /* Modify this basic chi-squared value as described in the Stutski &
385    Gusten paper. */
386       if( ndim == 1 ) {
387          t = ( cupidGC.beam_sq > 0.0 ) ? x0_off*x0_off/cupidGC.beam_sq : 0.0;
388       } else {
389          t = ( cupidGC.beam_sq > 0.0 ) ?
390                ( x0_off*x0_off + x1_off*x1_off )/cupidGC.beam_sq : 0.0;
391          if( ndim == 3 && cupidGC.velres_sq > 0.0 ) t += v_off*v_off/cupidGC.velres_sq;
392       }
393       chisq += cupidGC.sa*pdiff*pdiff + cupidGC.sc4*t + back_term;
394 
395 /* Store more diagnostic info */
396       if( cupidGC.nf == 1 ) {
397          pim = cupidGC.initmodel;
398          pm = cupidGC.model;
399          for( iel = 0; iel < cupidGC.nel; iel++ ) *(pim++) = *(pm++);
400       }
401       cupidGC.chisq = chisq;
402 
403    }
404 
405 /* Select or calculate the required return value.  If the chi squared
406    value itself is required, just return the value found above. */
407    if( what < 0 ) {
408       ret = chisq;
409 
410         cupidGCDumpF( MSG__DEBUG3, NULL, 0, NULL, NULL, status );
411 
412          msgSeti( "NF", cupidGC.nf );
413          msgOutif( MSG__DEBUG3, "", "   Fit attempt ^NF:", status );
414 
415          msgSetd( "C", ret );
416          msgOutif( MSG__DEBUG3, "", "      Chi-squared: ^C", status );
417 
418          msgSetd( "V", par[ 0 ] );
419          msgOutif( MSG__DEBUG3, "", "      Peak intensity: ^V", status );
420          msgSetd( "V", par[ 1 ] );
421          msgOutif( MSG__DEBUG3, "", "      Constant background: ^V", status );
422          msgSetd( "V", par[ 2 ] );
423          msgOutif( MSG__DEBUG3, "", "      Centre on 1st axis: ^V", status );
424          msgSetd( "V", par[ 3 ] );
425          msgOutif( MSG__DEBUG3, "", "      FWHM on 1st axis: ^V", status );
426 
427          if( ndim > 1 ) {
428             msgSetd( "V", par[ 4 ] );
429             msgOutif( MSG__DEBUG3, "", "      Centre on 2nd axis: ^V", status );
430             msgSetd( "V", par[ 5 ] );
431             msgOutif( MSG__DEBUG3, "", "      FWHM on 2nd axis: ^V", status );
432             msgSetd( "V", par[ 6 ] );
433             msgOutif( MSG__DEBUG3, "", "      Position angle: ^V", status );
434 
435             if( ndim > 2 ) {
436                msgSetd( "V", par[ 7 ] );
437                msgOutif( MSG__DEBUG3, "", "      Centre on vel axis: ^V", status );
438                msgSetd( "V", par[ 8 ] );
439                msgOutif( MSG__DEBUG3, "", "      FWHM on vel axis: ^V", status );
440                msgSetd( "V", par[ 9 ] );
441                msgOutif( MSG__DEBUG3, "", "      Vel gradient on 1st axis: ^V", status );
442                msgSetd( "V", par[ 10 ] );

      
   


class GaussClumps:

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
      # which do not influence the Gaussian model, status ). 
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
  
   def self.optimize():
      # Unpack used parameters
      maxnf=self.par['MAXNF']
      wwidth=self.par['WWIDTH']
      wmin=self.par['WMIN']
      rms=self.par['RMS']
      beam=self.par['BEAM']
      velres=self.par['VELORES']
      beamfwhm=self.par['FWHMBEAM']
      self.bfsq=beamfwhm*beamfwhm
      self.velsq=velres*velres

      # Gaussian Window

      # The factor which scales the FWHM on each axis to the half-width of the
      # section of the data array to be be fitted. 
      beta=0.5*wwidth*npsqrt(-np.log( wmin )/ np.log( 2.0 ) );
      lb,ub=self.res.index_from_window(self.cval,beta*self.fobs)
      # Store the data normalised to the
      # RMS noise level. Also calculate and store the Gaussian weight for the
      # pixel.
      y=self.res.get_slice(lb,ub).ravel()
      feat=self.res.get_features(lb,ub)
      wpos=np.array([self.cval[0],self.cval[1]])
      wstd=np.array([self.fobs[0],self.fobs[1]])*wwidth
      wfreq=self.cval[2]
      wfwhm=self.fobs[2]*wwidth
      (wmu,wP)=flx.clump_to_gauss(wpos,wstd,0,wfreq,wfwhm,np.array[0,0])
      # Normalise the weights to a maximum value of 1.0 and set to zero any weights
      # which are lower than the user supplied lower limit.
      we=flx.create_gauss(wmu,wP,feat,1.0)
      we[we < wmin]=0.0
      
      # Normalise all other data values in the guess structure and in the 
      # array to the RMS noise level.
      y=y/rms
      self.ymax /= rms;
      guess[1] /= rms;
      guess[0] /= rms;

      # Number of invocations of the function
      self.nf=0
      
      # TODO: Check and understand this
      # Get the factor by which to correct the peak amplitude of the model to
      # take account of the smoothing by the instrumental beam.
      t = guess[3]*guess[3]
      dx_sq = self.bfsq + t;
      peakfactor = t/dx_sq;
      t = guess[5]*guess[5];
      dx_sq = self.bfsq + t;
      peakfactor *= t/dx_sq;
      t = guess[8]*guess[8];
      dx_sq = self.velsq + t;
      peakfactor *= t/dx_sq;

      # Do the correction.
      if  peakfactor > 0.0:
         guess[0] /= sqrt(peakfactor);
      
      if self.fixback:
         np.delete(guess,[1])

      self.reuse=False
      # Optimize at last!
      res=fmin_bfgs(chi2, guess,fprime=jac_chi2, args=self)
      
      # TODO: Repeat when failed!
 

   def setInit(self,niter):
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
      (self.fobs[0],off[0]) = self.profWidth(self.res,self.imax,0) 
      (self.fobs[1],off[1]) = self.profWidth(self.res,self.imax,1) 
      (self.fobs[2],off[2]) = self.profWidth(self.res,self.imax,2) 
      fbeam=0.5*(self.fobs[0]  + self.fobs[1])/beamfwhm
      if fbeam < 1.0: 
         fbeam=1.2
      self.fobs[0] = fbeam*beamfwhm
      self.fobs[1] = fbeam*beamfwhm
           
      # Store the Guessed model 
      self.cval=self.res.index_to_wcs(self.imax)
      guess[2]=self.cval[0]
      guess[4]=self.cval[1]
      guess[7]=self.cval[2]
      # Find the initial guess at the intrinsic FWHM (i.e. the FWHM of the
      # clump before being blurred by the instrument beam). Do the same for 
      # the second axis. Assume zero rotation of the elliptical clump shape.
      guess[3]=np.sqrt(fbeam*fbeam- 1.0 )*beamfwhm
      guess[5]=guess[4]
      guess[6]=0.0
      # Now do the same for the third (velocity) axis if necessary. Assume
      # zero velocity gradient
      fvel=self.fobs[2]/velres
      guess[8]=np.sqrt(fvel*fvel- 1.0 )*velres
      guess[9]=0.0
      guess[10]=0.0
      
      # Store the mean of the background estimates, and the peak value. Noise
      # will result in the peak data value being larger than the peak clump value
      # by about the RMS noise. Therefore, reduce the peak value by the RMS.
      guess[1] = off.sum()/3
      guess[0] = self.ymax - guess[1] - rms

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
         guess[0] += guess[1];
         guess[1] = 0.0;
         self.fixback = True;
       self.guess=guess

   def fit(self,cube,verbose=False,use_meta=True):
      # Set the RMS, or automatically find an estimate for it
      if not self.par.has_key('RMS'):
         rms=cube.estimate_rms()
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
      self.res=cube.copy()

      # Initialise the number of clumps found so far.
      iclump = 0

      # Indicate that no peaks have been found below the lower threshold for clump
      # peak values, or below the lower area threshold. 
      peaks_below = 0
      area_below = 0
 
      # Initialise the variables used to keep track of the mean and standard
      # deviation of the most recent "npeak" fitted peak values. 
      mean_peak = 0.0;
      sigma_peak = 0.0; 
      # The value most recently added to "peaks"
      new_peak = 0.0;
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
      sumdata = res.get_flux()
      
      # peaks contains the last npeaks... 
      peaks=np.zeros(npeaks)
      
      # Loop round fitting a gaussian to the largest remaining peak in the
      # residuals array. */
      while iterate:
         # Report the iteration number to the user if required.
         niter+=1
         if verbose:
            log.info("Iteration: "+str(niter))
         # Find the cube index of the element with the largest value in the residuals cube.
         # imax: Index of element with largest residual
         (fmax,imax) = res.max()
 
         # Finish iterating if all the residuals are bad, or if too many iterations
         # have been performed since the last succesfully fitted clump. 
         if np.isnan(imax):
            iterate = False;
            niter-=1
            if verbose:
               log.info("There are no good pixels left to be fitted.")
               continue
         elif  nskip > maxskip:
            iterate = False;
            niter-=1
            if verbose:
               log.info("The previous",maxskip,"fits were unusable.")
               continue
         # If not, make an initial guess at the Gaussian clump parameters centred on the current peak.
         self.setInit()
         
         # Find the best fitting parameters, starting from the above initial guess.
         (found,clump,chisq)=self.optimize()
         # If no fit could be performed, then found = False
         if found:
            # Skip this fit if we have an estimate of the standard deviation of the
            # "npeaks" most recent clump peak values, and the peak value of the clump
            # just fitted is a long way (more than NSIGMA standard deviations) from the
            # peak value of the previously fitted clump. Also skip it if the peak
            # value is less than the "mlim" value.
            if (peaks.size == 0 or iclump < npeak or np.abs(clump[0] - new_peak) < nsig*sigma_peak ) and clump[0] > mlim]:

               # Record the new peak value for use with the next peak, and update the
               # standard deviation of the "npeak" most recent peaks. These values are
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
               nskip = 0;

               # Remove the model fit (excluding the background) from the residuals.
               # This also creates data values asociated with the clumps
               # The standard deviation of the new residuals is returned. */
               (result,area,sumclumps)=self.updateResults(clump,imax,mean_peak,sumclumps) # check check check
               #cupidGCUpdateArrays( type, res, ipd, el, ndim, dims,
               #                    x, rms, mlim, imax, peak_thresh, slbnd,
               #                    &ret, iclump, excols, mean_peak,
               #                    maxbad, &area, &sumclumps, status );

               # TODO: implement this!
               # Display the clump parameters on the screen if required. */
               #cupidGCListClump( iclump, ndim, x, chisq, slbnd, rms, status );

               # If this clump has a peak value which is below the threshold, increment
               # the count of consecutive clumps with peak value below the threshold.
               # Otherwise, reset this count to zero.
               if clump[0] < peak_thresh:
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
               res[imax]=np.nan;
               new_peak = 0.5*(new_peak + clump[0]);
               nskip++;
               if verbose:
                  log.info("Clump rejected due to aberrant peak value. Ignoring Pixel...")

         # Tell the user if no clump could be fitted around the current peak
         # pixel value 
         else:
            nskip++;
            if verbose:
               log.info("No clump fitted (optimization falied). Ignoring Pixel...")
            # Set the specified element of the residuals array bad if no fit was
            # performed. This prevents the any subsequent attempt to fit a Gaussian
            # to the same peak value.
            res[imax]=np.nan;
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
      return result

#def chi2(model,features,values,w,value_max,feature_max,res_vect,s_vect):
#   sys.stdout.write('.')
#   sys.stdout.flush()
#   su=values - gauss_eval(features,to_gauss(model))
#   nf=len(su) - 11;
#   t1 = (np.square(su)*w).sum()/nf
#   t2 = (np.exp(-su)).sum()/nf
#   t3 = np.square((model[2]-feature_max[0])/res_vect[0]) + np.square((model[3]-feature_max[1])/res_vect[1]) + np.square((model[4]-feature_max[2])/res_vect[2])
#   t4 = np.square(model[0]+model[1] - value_max)
#   return(t1 + s_vect[0]*t2 + s_vect[1]*t3 + s_vect[2]*t4)
#
#def jac_chi2(model,features,values,w,value_max,feature_max,res_vect,s_vect):
#   sys.stdout.write('*')
#   sys.stdout.flush()
#
#   # Unpack values
#   (a, b, x0, y0, v0, phi, sx, sy, sv, dvx, dvy) = model
#   (a,b,mu,L)=to_gauss(model)
#   fit=gauss_eval(features,(a,b,mu,L))
#   su=values - fit
#   nf=len(su) - 11;
#   # Compute basic term
#   basic_term = (-2*su*w + s_vect[0]*np.exp(-su))/nf
#   exp_term = (fit - b)/a
#   jaco=np.empty(11)
#   # partial derivate w.r.t. a
#   extra_sa = 2*s_vect[2]*(a + b - value_max)
#   jaco[0]=(basic_term*exp_term).sum()  + extra_sa
#   # partial derivate w.r.t. b
#   jaco[1]=basic_term.sum() + extra_sa
#   # partial derivate w.r.t. x0,y0 and v0
#   C=np.empty_like(features)
#   C[0]=features[0] - mu[0]
#   C[1]=features[1] - mu[1]
#   C[2]=features[2] - mu[2]
#   V=-L.dot(C) # Derivate of the quadratic form
#   extra_sc = 2*s_vect[1]*(mu - feature_max)/(res_vect*res_vect)
#   # compose x0
#   jaco[2]=(a*basic_term*exp_term*V[0]).sum() + extra_sc[0]
#   # compose y0
#   jaco[3]=(a*basic_term*exp_term*V[1]).sum() + extra_sc[1]
#   # compose v0
#   jaco[4]=(a*basic_term*exp_term*V[2]).sum() + extra_sc[2]
#   ## partial derivate w.r.t. phi,sx,sy,sv,dvx,dvy
#   sphi=np.sin(phi)
#   s2phi=np.sin(2*phi)
#   sphi2=np.square(sphi)
#   cphi=np.cos(phi)
#   c2phi=np.cos(2*phi)
#   cphi2=np.square(cphi)
#   sx2=np.square(sx)
#   sy2=np.square(sy)
#   sv2=np.square(sv)
#   sx3=np.power(sx,3)
#   sy3=np.power(sy,3)
#   sv3=np.power(sv,3)
#   dvx2=np.square(dvx)
#   dvy2=np.square(dvy)
#   D_phi=(1.0/sy2 - 1.0/sx2)*np.array([[2*sphi*cphi,c2phi,0],[c2phi,-2*sphi*cphi,0],[0,0,0]])
#   V=-(C*(D_phi.dot(C))).sum(axis=0)/2.0
#   jaco[5]=(a*basic_term*exp_term*V).sum()
#   D_sx=(1.0/sx3)*np.array([[-2*cphi2,s2phi,0],[s2phi,-2*sphi2,0],[0,0,0]])
#   V=-(C*(D_sx.dot(C))).sum(axis=0)/2.0
#   jaco[6]=(a*basic_term*exp_term*V).sum()
#   D_sy=(1.0/sy3)*np.array([[-2*sphi2,-s2phi,0],[-s2phi,-2*cphi2,0],[0,0,0]])
#   V=-(C*(D_sy.dot(C))).sum(axis=0)/2.0
#   jaco[7]=(a*basic_term*exp_term*V).sum()
#   D_sv=(2.0/sv3)*np.array([[-dvx2,-dvx*dvy,dvx],[-dvx*dvy,-dvy2,dvy],[dvx,dvy,-1]])
#   V=-(C*(D_sv.dot(C))).sum(axis=0)/2.0
#   jaco[8]=(a*basic_term*exp_term*V).sum()
#   D_dvx=(1.0/sv2)*np.array([[2*dvx,dvy,-1],[dvy,0,0],[-1,0,0]])
#   V=-(C*(D_dvx.dot(C))).sum(axis=0)/2.0
#   jaco[9]=(a*basic_term*exp_term*V).sum()
#   D_dvy=(1.0/sv2)*np.array([[0,dvx,0],[dvx,2*dvy,-1],[0,-1,0]])
#   V=-(C*(D_dvy.dot(C))).sum(axis=0)/2.0
#   jaco[10]=(a*basic_term*exp_term*V).sum()
#   return jaco
#
#def next_clump(cube,syn,params):
#   # Non-blocking plot 
#   plt.ion() 
#   plt.clf() 
#   
#   (value_max,feature_max) = cube.max()
#   b=params['rms'] # need to be a local minima... not good
#   a=value_max - b
#   
#   # Initial guess: position and orientation
#   x0=feature_max[0]
#   y0=feature_max[1]
#   v0=feature_max[2]
#   phi=0
#   
#   # Initial guess: variances
#   res_vect=np.array([params['beam_size'],params['beam_size'],params['spe_res']])
#   s_vect=np.array([params['s0'],params['sc'],params['sa']])
# 
#   # Initial guess: redshift
#   dvalp=0
#   dvdel=0
#  
#   # Compute the weight vector and the feature space
#   w_sigmas=params['weight_deltas']*res_vect # several times the resolution
#   (features,sc_index) = cube.feature_space(feature_max,2*w_sigmas) 
#   w_shape=(1,0,x0,y0,v0,0,w_sigmas[0],w_sigmas[1],w_sigmas[2],0,0)
#   w=gauss_eval(features,to_gauss(w_shape),False)
# 
#   #Plot current cube
#   plt.subplot(2, 3, 4)
#   plt.imshow(cube.stack())
#   rect=plt.Rectangle((sc_index[0],sc_index[2]),sc_index[1]-sc_index[0]+1,sc_index[3]-sc_index[2]+1,alpha=1, facecolor='none')
#   plt.gca().add_patch(rect)
#   
#   #Plot current subcube
#   plt.subplot(2, 3, 1)
#   plt.imshow(cube.stack(sc_index))
#
#   
#   # Ravel the values
#   values=cube.ravel(sc_index)
#   
#   idx=(values-a/2).argmin()
#   sq2l2=np.sqrt(2*np.log(2))
#   sx=np.abs(features[0][idx]-feature_max[0])/sq2l2
#   sy=np.abs(features[1][idx]-feature_max[1])/sq2l2
#   sv=np.abs(features[2][idx]-feature_max[2])/sq2l2
#   # Compile first guess
#   guess=np.array([a,b,x0,y0,v0,phi,sx,sy,sv,dvalp,dvdel])
#   print "GUESS"
#   print guess
#   
#   # Pack all args of chi2
#   chi2_args=(features,values,w,value_max,feature_max,res_vect,s_vect)
#   
#
#   print "grad error", check_grad(chi2,jac_chi2,guess,features,values,w,value_max,feature_max,res_vect,s_vect)
#   #print "PREV= ", chi2_args
#   # Standarize everything
#   #(std_args,tr_vect)=standarize(chi2_args)
#   #print "STAND= ", std_args
#   #print "tr= ", tr_vect
#   #s_guess=guess/tr_vect
#   # OPTIMIZE
#   #res = minimize(chi2,guess,jac=jac_chi2,method='CG',args=chi2_args)
#   res=fmin_bfgs(chi2, guess,fprime=jac_chi2, args=chi2_args)
#   #res=fmin_bfgs(chi2, guess,args=chi2_args)
#   #res = minimize(chi2,guess,jac=jac_chi2,method='BFGS',args=chi2_args,tol=1e-30)
#   print
#   print "res =", res
#   print
#   print "clump =", res  #*tr_vect
#   sys.stdout.flush()
#   # Clump values
#   #print "AND NOW THE GAUSS"
#   #print to_gauss(res.x) 
#   val_fit=gauss_eval(features,to_gauss(res),False)
#   fit_cube=cube_data_unravel(val_fit,sc_index)
#   # Remove clump from the real cube 
#   cube.add(-fit_cube,sc_index)
#   # Add clump to the synthetic cube
#   syn.add(fit_cube,sc_index)
#   
#   #Plot current cube
#   plt.subplot(2, 3, 5)
#   plt.imshow(cube.stack())
#   plt.gca().add_patch(rect)
#   
#   #Plot current subcube
#   plt.subplot(2, 3, 2)
#   plt.imshow(cube.stack(sc_index))
#   
#   #Plot clump
#   plt.subplot(2, 3, 3)
#   plt.imshow(cube_data_stack(fit_cube))
#
#   #Plot synthetic
#   plt.subplot(2, 3, 6)
#   plt.imshow(syn.stack())
#   plt.gca().add_patch(rect)
#
#   #END
#   #M=clu.reshape((v_ub-v_lb+1,y_ub-y_lb+1,x_ub-x_lb+1))
#   
#   #cube.data[v_lb:v_ub+1,y_lb:y_ub+1,x_lb:x_ub+1] -= M
#   #syn[v_lb:v_ub+1,y_lb:y_ub+1,x_lb:x_ub+1] += M
#   
#   # Matrices for displaying results (NOT ALGORITHMIC CODE)
#   #ma=cube.data[v_lb:v_ub+1,y_lb:y_ub+1,x_lb:x_ub+1].sum(axis=0)
#   #spe=cube.data[v_lb:v_ub+1,y_lb:y_ub+1,x_lb:x_ub+1].sum(axis=(1,2))
#   #prof=M.sum(axis=0)
#   #vmin=ma.min()
#   #vmax=ma.max()
#   #plt.clf() 
#   #plt.subplot(2, 3, 1)
#   #plt.imshow(ma,vmin=vmin,vmax=vmax)
#   #plt.subplot(2, 3, 3)
#   #plt.imshow(prof)
#   #plt.subplot(2, 3, 2)
#
#  
#
#   # SHOW results (NOT ALGORITHMIC CODE)
#   #spe2=cube.data[v_lb:v_ub+1,y_lb:y_ub+1,x_lb:x_ub+1].sum(axis=(1,2))
#   #plt.imshow(cube.data[v_lb:v_ub+1,y_lb:y_ub+1,x_lb:x_ub+1].sum(axis=0),vmin=vmin,vmax=vmax)
#   #plt.subplot(2, 3, 6)
#   #plt.imshow(syn.sum(axis=0))
#   #plt.subplot(2, 3, 5)
#   #plt.imshow(cube.data.sum(axis=0))
#   #plt.gca().add_patch(plt.Rectangle((x_lb,y_lb),x_ub-x_lb+1,y_ub-y_lb+1,alpha=1, facecolor='none'))
#   #plt.subplot(2, 3, 4)
#   #plt.plot(spe,"b")
#   #plt.plot(spe2,"r")
#   plt.show() 
#   plt.pause(0.01)
#
#   # Return the clump parameters
#   return res
#
#
#
#def compute_rms(data):
#   res=data[data < 0]
#   fin=(res*res).sum()/len(res)
#   return np.sqrt(fin)
#
#def gauss_clumps(orig_cube,params):
#   cube=copy.deepcopy(orig_cube)
#   syn=copy.copy(orig_cube)
#   syn.data=np.empty_like(cube.data)
#   C=[]
#   stop=False
#   norm=cube.data.mean()
#   print "Initial Sum", norm
#   params['beam_size']=abs(float(cube.meta['BMIN']))
#   params['spe_res']=abs(float(cube.meta['CDELT3']))
#   params['rms']=compute_rms(cube.data)
#   print params['rms']
#   while not stop:
#      theta=next_clump(cube,syn,params)
#      C.append(theta)
#      norm=cube.data.mean()
#      print "Sum", norm
#      stop = (norm < params['threshold'])
#   
#   # SHOW RESULTS
#   plt.clf() 
#   plt.subplot(1, 3, 2)
#   plt.imshow(cube.stack())
#   plt.subplot(1, 3, 1)
#   plt.imshow(orig_cube.stack())
#   plt.subplot(1, 3, 3)
#   plt.imshow(syn.stack())
#   plt.show()
#   plt.pause(100)
#   print "RESULTS"
#   print C
#   return C
