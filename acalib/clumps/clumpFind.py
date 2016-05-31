import copy
import sys
import ca
import numpy as np
import matplotlib.pyplot as plt


class FellWalker:

   def __init__(self):
      self.defaultParams()

   def defaultParams(self):
      self.par = dict()

      """ Generic parameters  """
      # Spectral resoluion in pixels
      self.par['VELORES'] = 2.0
      # Beam resolution in pixels
      self.par['FWHMBEAM'] = 2.0
      # Maximum Clumps
      self.par['MAXCLUMPS'] = sys.maxint
      # The lower threshold for clump values to a user-specified multiple of the RMS noise.
      self.par['THRESH'] = 1.0

      """ Specific ClumpFind parameters """
      #Numbers of axis to consider to compute the neighboring pixels
      self.par['NAXIS'] = 3


   def create_caa(self, data):
      caa = np.zeros_like(data.data).astype(np.int)
      #Check invalid pixels (below threshold)
      rms = self.par["RMS"]
      threshold = self.par['THRESH']*rms
      #Aditionally, NaN valued pixels are set as unusable -> filled(1)
      mask = np.array((data<threshold).filled(1))
      caa[mask] = -1
      return caa


   def compute_levels(self, maxv, minv, rms):
      #Get the lowest contour level using twice the RMS as the default.
      clow = 2*self.par['RMS']

      #Report an error if the lowest contour level is below the minimum value
      #in the data array.
      if clow < minv and verbose:
         log.warning('The supplied lowest contour level is below the minimum value in the data')

      #Report an error if the lowest contour level is above the maximum value
      #in the data array
      if clow >= maxv and verbose:
         log.warning('The supplied lowest contour level is above the maximum value in the data')

      #Use 2*RMS as the contour delta increment.
      cdelta = 2*self.par['RMS']

      #Find the number of levels needed for this delta
      nlevels = int((maxv-clow)/cdelta)+1
      self.par['NLEVELS'] = nlevels

      #Allocate the array and fill it with the appropriate contour levels
      ret = np.arange(clow, maxv, cdelta)
      ret = ret[::-1]
      return ret

   def scan(self, data, caa, clevel, naxis):
      return


   # Back compatibility function
   def fit(self, orig_cube, verbose=False):
       log.warning("The .fit() interface will be removed soon, please use .run()")
       return self.run(cube, verbose)

   def run(self, orig_cube, verbose=False):
      cube = orig_cube.copy()
      syn = orig_cube.empty_like()

      # Set the RMS, or automatically find an estimate for it
      if not self.par.has_key('RMS'):
         rms = cube.estimate_rms()
         self.par['RMS'] = rms

      #Number of axis to compute neighborhood
      naxis = self.par['NAXIS']

      #cube.data is masked array
      data = cube.data

      """
      Fill the supplied caa array with -1 for all pixels which are below the
      threshold, or are bad. Fill all other pixels with zero to indicate that
      the pixel is "usable but not yet checked".
      """
      print "Creating CAA"
      caa=self.create_caa(data)

   
      #Initialize dictionary which describe clumps
      self.clump = dict()

      #Find the largest and mallest good data values in the supplied array.
      maxv = np.max(data)
      minv = np.min(data)

      #Get the contour levels at which to check for clumps
      levels = self.compute_levels(maxv, minv, rms)

      #Initialise the largest data value in the remaining unassigned pixels.
      self.maxrem = maxv

      #Loop round all contour levels.
      for clevel in levels:
         """
         Scan the data array at a new contour level. This extends clumps found
         at a higher contour level, and adds any new clumps found at this contour
         level. New clumps are stored at the end of the returned array. If the
         current contour level is higher than the maximum of the remaining
         unassigned pixel values, there is no point in doing this scan since it
         will find no pixels.
         """
         if clevel < self.maxrem:
            clumps = self.scan(data, caa, clevel, naxis)

         elif:
            log.warning('No pixels found at this contour level.')

   #return caa and clump structures
   return caa,clump
