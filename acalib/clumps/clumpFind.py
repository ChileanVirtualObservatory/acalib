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

      #NaN valued pixels are set as unusable -> filled(1)
      mask = np.array(data.filled(1))
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

   def neighborhood(self, pos, caa, naxis):
      #dimensions of caa array
      dims = caa.shape
      #actual position
      (x,y,z) = pos
      #clump labels of neighbors
      ret = list()

      #TODO
      if naxis==1:
         return ret

      #TODO
      elif naxis==2:
         return ret

      elif naxis==3:
         for i in [-1,0,1]:
            for j in [-1,0,1]:
               for k in [-1,0,1]:
                  if i==j==k==0: continue
                  neigh = (x+i,y+j,z+k)
                  #out condition 1
                  out1 = neigh[0]<0 or neigh[1]<0 or neigh[2]<0
                  #out condition 2
                  out2 = neigh[0]>=dims[0] or neigh[1]>=dims[1] or neigh[2]>=dims[2]

                  if out1 or out2: continue
                  elif caa[neigh]<=0: continue
                  else: ret.append(caa[neigh])
      return ret.sort()

   def scan(self, data, caa, clevel, naxis):
      """
      Note the next available index value. This value is saved now so that
      we can differentiate later between PixelSets created at the current
      contour level and those created at higher contour levels (PixelSets
      created at higher contour levels will have indices less than hindex).
      """
      hindex = self.index

      """
      Scan the data array for good pixels which are above (or at) the supplied
      contour level and have not yet been assigned to a PixelSet (i.e. have a
      null index in the caa array). Keep a check on whether the pixel is
      an edge pixel or not (TODO).
      """
      mask = np.logical_and(data>=clevel, caa==0)
      positions = np.array(np.where(mask)).T

      for pos in positions:
         """
         When such a pixel is found, have a look at its immediate neighbours and
         see if any of them have already been assigned to a PixelSet. If they
         have, identify which PixelSet they are assigned to. We distinguish two
         different types of PixelSets; those which were identified at this contour
         level, and those which were identified at higher contour levels.
         """
         neighbors = self.neighborhood(pos, caa, naxis, hindex)



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

      #Initialise the index used to identify the next contiguous set of
      #pixels found
      self.index = 1

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
