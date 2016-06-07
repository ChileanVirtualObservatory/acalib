import copy
import sys
import ca
import numpy as np
import matplotlib.pyplot as plt
from math import sqrt

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
      self.par['NAXIS'] = 2


   def create_caa(self, data):
      caa = np.zeros_like(data.data).astype(np.int)
      #NaN valued pixels are set as unusable -> filled(1)
      mask = np.array(data.filled(1))
      caa[mask] = -1
      return caa

   #square distance between p0 and p1
   def dist(self,  p0, p1):
      if len(p0)==len(p1)==2:
         return (p0[0]-p1[0])**2 + (p0[1]-p1[1])
      elif len(p0)==len(p1)==3:
         return (p0[0]-p1[0])**2 + (p0[1]-p1[1])**2 + (p0[2]-p1[2])**2
      else: return None


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

   def neighborhood(self, pos, caa, naxis, hindex):
      #dimensions of caa array
      dims = caa.shape
      ndim = len(dims)

      #number of PixelSets defined at the current contour 
      #level which adjoin the pixel specified by "pos"
      n1 = 0

      #The first "n1" elements of this list will be returned holding
      #the indices of all PixelSets defined at the current contour level
      #which adjoin the pixel specified by "pos".
      i1 = list()

      #The lowest index of any PixelSets defined at the current contour
      #level which adjoin the pixel specified by "pos"
      i11 = float('inf')

      #The number of PixelSets defined at a higher contour level which 
      #adjoin the pixel specified by "pos"
      n2 = 0

      #The index of a PixelSet defined at a higher contour level which 
      #adjoins the pixel specified by "pos". If there is more than one 
      #such PixelSet, the returned PixelSet is the one which has the closest 
      #peak to the pixel being tested.
      i12 = None


      #2D data case
      if ndim==2:
         #actual position
         (x,y) = pos
         #clump labels of neighbors
         ret = list()

         #variating one dimension
         if naxis=1:
            for i in [-1,1]:
               neigh = (x+i,y)
            for j in [-1,1]:
               neigh = (x,y+j)

         #variating two dimensions
         elif naxis=2:
            for i in [-1,0,1]:
               for j in [-1,0,1]:
                  if i==j==0: continue
                  neigh = (x+i,y+j)
                  nindex = caa[neigh]
                  #out condition 1
                  out1 = neigh[0]<0 or neigh[1]<0
                  #out condition 2
                  out2 = neigh[0]>=dims[0] or neigh[1]>=dims[1]

                  if out1 or out2: continue
                  elif nindex <= 0: continue
                  else:
                     if nindex>=hindex:
                        n1+=1
                        i1.append(nindex)
                        if nindex < i11:
                           i11 = nindex
                     elif nindex<hindex:
                        n2+=1

      #3D data case
      if ndim==3:
         #actual position
         (x,y,z) = pos
         #clump labels of neighbors
         ret = list()

         #TODO
         if naxis==1:
            for i in [-1,1]:
               neigh = (x+i,y,z)
            for j in [-1,1]:
               neigh = (x,y+j,z)
            for k in [-1,1]:
               neigh = (x+y+z+k)
            return ret

         #TODO
         elif naxis==2:
            return ret

         #variating tree dimensions
         elif naxis==3:
            for i in [-1,0,1]:
               for j in [-1,0,1]:
                  for k in [-1,0,1]:
                     if i==j==k==0: continue
                     neigh = (x+i,y+j,z+k)
                     nindex = caa[neigh]
                     #out condition 1
                     out1 = neigh[0]<0 or neigh[1]<0 or neigh[2]<0
                     #out condition 2
                     out2 = neigh[0]>=dims[0] or neigh[1]>=dims[1] or neigh[2]>=dims[2]

                     if out1 or out2: continue
                     elif nindex<=0: continue
                     else:
                        if nindex>=hindex:
                           n1 += 1
                           i1.append(nindex)
                           if nindex < i11:
                              i11 = nindex

                        elif nindex<hindex:
                           n2 += 1
      return (n1, i11, i1, n2, i12)

   def scan(self, data, caa, clevel, naxis):
      """
      Note the next available index value. This value is saved now so that
      we can differentiate later between PixelSets created at the current
      contour level and those created at higher contour levels (PixelSets
      created at higher contour levels will have indices less than hindex).
      """
      hindex = self.index

      """
      Structure containing the indexes of clumps identified at higher contour
      levels, that adjoins each of the clumps identified at the current level
      """
      adjoins = dict()

      """
      Scan the data array for good pixels which are above (or at) the supplied
      contour level and have not yet been assigned to a PixelSet (i.e. have a
      null index in the caa array). Keep a check on whether the pixel is
      an edge pixel or not (TODO).
      """
      mask = np.logical_and(data>=clevel, caa==0)
      positions = np.array(np.where(mask)).T
      positions = map(tuple, positions)

      for pos in positions:
         """
         When such a pixel is found, have a look at its immediate neighbours and
         see if any of them have already been assigned to a PixelSet. If they
         have, identify which PixelSet they are assigned to. We distinguish two
         different types of PixelSets; those which were identified at this contour
         level, and those which were identified at higher contour levels.
         """
         (n1, i11, i1, n2, i12) = self.neighborhood(pos, caa, naxis, hindex)

         """
         If none of the neighbours of this pixel are assigned to a PixelSet which
         was identified at this contour level, then we start a new PixelSet.
         """
         if n1==0:
            self.clumps[self.index] = [pos]
            self.peaks[self.index] = pos
            caa[pos] = self.index
            self.index += 1

         """
         If one or more of the neighbours of this pixel are assigned to PixelSets
         which were identified at this contour level, then add this pixel into
         the PixelSet with the lowest index
         """
         else:
            (self.clumps[i11]).append(pos)
            if data[pos]>self.peaks[i11]: self.peaks[i11] = pos
            caa[pos] = i11 
            """
            If this pixel touches other PixelSets identified at this contour level,
            then transfer the pixels contained in them all into the PixelSet with
            lowest index
            """
            if n1>1:
               #removing repeated indexes and i11 index
               i1 = set(i1)
               i1.remove(i11)

               #updating clumps dict, peaks dict and caa 
               for ind in i1:
                  positions = self.clumps.pop(ind)
                  peak_pos = self.peaks.pop(ind)
                  self.clumps[i11] += positions
                  if data[peak_pos]>self.peaks[i11]: self.peaks[i11] = peak_pos 
                  for pos in positions:
                     caa[pos] = i11

      """
      Now check each of the new PixelSets created above. Ordering
      indexes of new clumps, in a sequential way
      """
      if self.index > hindex:
         seq_ind = hindex
         for clumps.keys()

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

   
      """
      Initialize structures wich describe clumps:
      clumps : {clump_index:[list of positions]}
      peaks : {clump_index: position of peak}
      """
      self.clumps = dict()
      self.peaks = dict()

      #Find the largest and smallest good data values in the supplied array.
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
            self.scan(data, caa, clevel, naxis)

         elif:
            log.warning('No pixels found at this contour level.')

   #return caa and clump structures
   return caa,clumps
