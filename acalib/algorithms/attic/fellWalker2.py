import copy
import sys
import ca
import numpy as np
import matplotlib.pyplot as plt
from astropy import log
import astropy.units as u
from astropy.nddata import *

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

      """ Specific FellWalker parameters """
      #Longest jump to a higher neighbouring pixel
      self.par['MAXJUMP'] = 6
      #The minimum dip parameter is specified as a multiple of the noise level in the data
      self.par['MINDIP'] = 3
      #Minimum size (in pixels) of clumps
      self.par['MINSIZE'] = 20
      #Lowest gradient which marks the start of the walk
      self.par['FLATSLOPE'] = 1.
      #Data values less than this are at "sea level"
      self.par['SEALEVEL'] = 1.

      """ cellulars automaton parameters """
      #Min fraction of neighbouring good pixels
      self.par['FRAC'] = 0.1
      #The value used to represent "on" pixels in input and output arrays
      self.par['ON'] = 0
      #The value used to represent "off" pixels in the output array (any
      #value not equal to "on" is treated as off in the input array).
      self.par['OFF'] = -1
      #If non-zero, then no output pixel will be set on if the
      #corresponding input pixel is not on.
      self.par['CENTRE'] = 1
      #Number of iterations of second celullar automaton.
      self.par['CLEANITER'] = 1


   def create_caa(self, data):
      caa = np.zeros(data.shape, dtype=np.int)
      #Check invalid pixels (below threshold)
      threshold = self.par['THRESH']*self.par["RMS"]
      #Aditionally, NaN valued pixels are set as unusable -> filled(1)
      if type(data)==np.ma.MaskedArray:
         mask = np.array((data<threshold).filled(1))
      else:
         mask = data<threshold
      caa[mask] = -1
      return caa


   def max_gradient(self, pos, data, caa):
      """
      We will now examine the 3x3x3 cube of neighbouring pixels to determine
      where we step to next. We will choose the neighbour with the greatest
      gradient.
      """
      max_pos=pos
      max_grad=0.
      dist=[np.sqrt(1),np.sqrt(2),np.sqrt(3)]
      shape=data.shape

      for i in range(-1,2):
         for j in range(-1,2):
            for k in range(-1,2):
               neigh=(pos[0]+i,pos[1]+j,pos[2]+k) #position of neighbour point
               out1=neigh[0]<0 or neigh[1]<0 or neigh[2]<0 #out condition 1
               out2=neigh[0]>=shape[0] or neigh[1]>=shape[1] or neigh[2]>=shape[2] #out condition 2
               d=dist[abs(i)+abs(j)+abs(k)-1] #distance to occupy

               if i==j==k==0:
                  #dont check it again
                  continue
               elif out1 or out2:
                  #don't go outside of cube
                  continue
               elif caa[neigh]==-1:
                  #don't check unusable pixels
                  continue
               else:
                  grad=(data[neigh]-data[pos])/d #calculate gradient
                  if grad>max_grad:
                     max_pos=neigh
                     max_grad=grad
      return max_pos

   def verify_peak(self, pos, data, caa):
      max_pos=pos
      max_val=data[pos]
      shape=data.shape
      maxJump=self.par['MAXJUMP']

      for i in range(-maxJump,maxJump+1):
         for j in range(-maxJump,maxJump+1):
            for k in range(-maxJump,maxJump+1):
               neigh=(pos[0]+i,pos[1]+j,pos[2]+k) #position of neighbour point
               out1=neigh[0]<0 or neigh[1]<0 or neigh[2]<0 #out condition 1
               out2=neigh[0]>=shape[0] or neigh[1]>=shape[1] or neigh[2]>=shape[2] #out condition 2

               if abs(i)<=1 and abs(j)<=1 and abs(k)<=1:
                  #don't check it again
                  continue
               elif out1 or out2:
                  #don't go outside of cube
                  continue
               elif caa[neigh]==-1:
                  #don't check unusable pixels
                  continue
               elif data[neigh]>max_val:
                  max_pos=neigh
                  max_val=data[neigh]
      return max_pos

   """
   return dictionaries peaks and cols of the form:
   peaks={clumpId:peak,...}
   cols={clumpId:{neighId:col,...},...}
   """
   def clump_structs(self, clump, data, caa):
      #Initialize peaks ans cols dictionaries
      peaks=dict()
      cols=dict()
      for clumpId in clump:
         cols[clumpId]=dict()
      shape=data.shape

      for clumpId,pixels in clump.items():
         peaks[clumpId]=0. #max value in current clump
         for pos in pixels:
            #First, find the peak of clump
            if data[pos]>peaks[clumpId]:
               peaks[clumpId]=data[pos]
            #Second, verify if it is a border clump pixel
            isBorder=False
            for i in range(-1,2):
               for j in range(-1,2):
                  for k in range(-1,2):
                     neigh=(pos[0]+i,pos[1]+j,pos[2]+k) #position of neighbour point
                     out1=neigh[0]<0 or neigh[1]<0 or neigh[2]<0 #out condition 1
                     out2=neigh[0]>=shape[0] or neigh[1]>=shape[1] or neigh[2]>=shape[2] #out condition 2
                     if i==j==k==0:
                        #don't check itself
                        continue
                     elif out1 or out2:
                        #don't go outside of cube
                        continue
                     elif caa[neigh]>0 and caa[neigh]!=clumpId:
                        #border clump pixel found!
                        isBorder=True
                        neighId=caa[neigh]

                        #Update for current clump (clumpId)
                        if neighId not in cols[clumpId]:
                           cols[clumpId][neighId]=data[pos]
                        elif data[pos]>cols[clumpId][neighId]:
                           cols[clumpId][neighId]=data[pos]

                        #Update for neighour (neighId)
                        if clumpId not in cols[neighId]:
                           cols[neighId][clumpId]=data[neigh]
                        elif data[neigh]>cols[neighId][clumpId]:
                           cols[neighId][clumpId]=data[neigh]
                        break
                  if isBorder: break
               if isBorder: break
      return peaks,cols


   def merge(self, clump, peaks, cols, caa, minDip):
      """
      Enter an iterative loop in which we join clumps that have a small dip
      between them. Finish when an interation fails to merge any clumps.
      """
      merge=True
      while merge:
         #Don't iterate again unless some merge occur
         merge=False
         for clumpId in clump:
            """
            Check all the clumps that adjoin clump "clumpId", and find the one
            that is separated from "clumpId" by the smallest dip (i.e. has the
            highest "col" ).
            """
            topcol=0.
            neighId=-1
            for tmpId,col in cols[clumpId].items():
               if col>topcol:
                  topcol=col
                  neighId=tmpId

            """
            If the peak value in clumpId is less than "mindip" above the col
            between clumpId and neighId, then we merge the clumps. */
            """
            if neighId!=-1 and peaks[clumpId]<topcol+minDip:
               #Merge it!
               merge=True
               break #dictionaries can't change size while iterating it
         if merge:
            """
            Tell the neighbours of clumpId that they are now neighbours of neighId
            instead. But only do this if they are not already neighbours of neighId.
            """
            #merging cols
            tmpCols=cols.pop(clumpId)
            del tmpCols[neighId]
            tmpCols.update(cols[neighId])
            del tmpCols[clumpId]
            cols[neighId]=tmpCols
            for tmpId in cols:
               if clumpId in cols[tmpId]:
                  cols[tmpId][neighId]=cols[tmpId].pop(clumpId)

            #delete clumpId from peaks and update new peak of merged clumps
            tmpPeak=peaks.pop(clumpId)
            peaks[neighId]=max(tmpPeak,peaks[neighId])

            #Update caa
            for pos in clump[clumpId]:
               caa[pos]=neighId

            #Merge pixels in clump dictionary
            clump[neighId]+=(clump.pop(clumpId))
      return clump,peaks,cols,caa


   def walkup(self, pos, path, pathv, data, caa):
      next_pos=self.max_gradient(pos,data,caa)

      if caa[next_pos]>=1:
         """
         If the next pixel on the route is already assigned to a
         clump, we now know what peak we are heading towards so use that clump
         index for the route so far, and abandon the rest of the walk.
         """
         path.append(next_pos)
         pathv.append(data[next_pos])
         return path,pathv

      elif next_pos!=pos:
         """
         Otherwise, move to the next pixel on the route.
         """
         path.append(next_pos)
         pathv.append(data[next_pos])
         return self.walkup(next_pos,path,pathv,data,caa)

      else:
         """
         If there is no upward route from the current pixel, (i.e. if it is a
         local maximum), we check the slightly more extended neighbourhood for
         a higher pixel. We do this because a local maximum could be just a
         noise spike.
         """
         local_max=next_pos
         new_max=self.verify_peak(local_max,data,caa)
         if local_max==new_max:
            """
            If the central pixel is the highest pixel in the more extended
            neighbourhood, we have reached a peak.
            """
            return path,pathv
         else:
            #Just a noise peak, keep walking up
            return self.walkup(new_max,path,pathv,data,caa)

         """
         The walk has now finished because we have either reached a peak which
         is the highest point in its neighbourhood, or we have joined an earlier
         route and thus know which peak we would get to if we were to carry on.
         """

   def verify_flat(self, path, pathv, caa, flatSlope, seaLevel):
      """
      Find the average gradient, if possible, over the last 4 steps, assuming
      each step is the same length. If it is above the user-supplied limit,
      indicate we have got the length of the initial flat section of the
      route.
      """
      valid=-1 #valid flag, to indicate from what pixel path is valid


      #In case of small paths (not in cupid version)
      if 1<len(path)<4:
            if pathv[0]>=seaLevel:
               valid=0
            else:
               avg=(pathv[-1]-pathv[0])/(len(path)-1)
               if avg>=flatSlope:
                  valid=0
      #In case of more than 4 steps in path
      else:
         for i in range(len(path)-3):
            if pathv[i]>=seaLevel:
               #valid from i-pixel
               valid=i
               break
            else:
               avg=(pathv[i+3]-pathv[i])/3
               if avg>=flatSlope:
                  #valid from i-pixel
                  valid=i
                  break

      #update path and flat
      if valid==-1:
         flat=path
         flatv=pathv
         path=list()
         pathv=list()
      else:
         flat=path[0:valid]
         flatv=pathv[0:valid]
         path=path[valid::]
         pathv=pathv[valid::]

      return path,pathv,flat,flatv

   # Back compatibility function
   #def fit(self, cube, verbose=False):
   #    log.warning("The .fit() interface will be removed soon, please use .run()")
   #    return self.run(cube, verbose=verbose)

   def run(self, data, rms, verbose=False):
      #cube = cube.copy()

      # Set the RMS, or automatically find an estimate for it
      #rms = estimate_rms(cube)
      self.par['RMS'] = rms

      #cube.data is masked array
      data = data.copy()

      """
      Fill the supplied caa array with -1 for all pixels which are below the
      threshold, or are bad. Fill all other pixels with zero to indicate that
      the pixel is "usable but not yet checked".
      """
      caa=self.create_caa(data)
      if verbose: log.info("CAA initialized successfully")

      """
      Use a cellular automata to remove small isolated groups of usable
      pixels from the above "ipa" array. This allocates new memory for the
      cleaned up array.
      """
      if verbose: log.info("[Stage 1] Removing small isolated regions")
      frac = self.par['FRAC']
      on = self.par['ON']
      off = self.par['OFF']
      centre = self.par['CENTRE']
      caa = ca.remove_isolate(caa, frac, on, off, centre)

      #Some constants
      noise = self.par['THRESH']*rms
      flatSlope = self.par['FLATSLOPE']*rms
      seaLevel = self.par['SEALEVEL']*rms


      clump=dict()
      shape=data.shape
      top_id=0 #top clump id

      """
      Scan through the caa array, looking for usable pixels which have not
      yet been assigned to a clump (i.e. have a value of zero in caa).
      """
      if verbose: log.info("[Stage 2] Walking up!")
      for i in range(shape[0]):
         for j in range(shape[1]):
            for k in range(shape[2]):
               # If unusable or already assigned, skip it
               if caa[i,j,k]!=0:
                  continue
               path=list() # Ascent path pixels positions
               pathv=list() #Ascent path pixel values
               pos=(i,j,k)
               path.append(pos)
               pathv.append(data[pos])


               """
               We now start walking away from this pixel up hill (up the steepest
               gradient). We store the vector index of all pixels visited on this walk
               in the route array. Whilst walking, we initially keep a record of the
               average gradient over the last 4 steps. We count how many steps we travel
               before we encounter an average gradient above the gradient limit.
               However, this initial flat section is only ignored if the walk starts from
               "sea level". Walks which start at a higher altitude are used in their entirety.
               """
               path,pathv=self.walkup(pos,path,pathv,data,caa)
               path,pathv,flat,flatv=self.verify_flat(path,pathv,caa,flatSlope,seaLevel)

               """
               We now assign the clump index found above to all the pixels visited on
               the route, except for any low gradient section at the start of the route,
               which is set unusable (-1). We ignore walks that were entirely on the
               coastal plain (indicated by a value of -2 for clump_index).
               """

               #if not empty
               if flat:
                  #set pixels as unusable
                  for pos in flat:
                     caa[pos]=-1

               #if after that, path is empty then continue
               if not path:
                  continue

               """
               If this peak has not already been assigned to a clump, increment the number
               of clumps and assign it the new clump index.
               """
               if caa[path[-1]]>0:
                  #Ascent path reach an existing path
                  path_id=caa[path[-1]]
                  for pos in path:
                     if pos==path[-1]:
                        #Already asigned to clump
                        continue
                     caa[pos]=path_id
                     clump[path_id].append(pos)
               else:
                  #A new ascent path
                  top_id+=1
                  path_id=top_id
                  clump[path_id]=list()
                  for pos in path:
                     caa[pos]=path_id
                     clump[path_id].append(pos)
      if verbose: log.info('Number of clumps found: {0}'.format(len(clump)))


      ###refine 1, removing small clumps
      if verbose: log.info("[Stage 3] Removing small clumps")
      minSize = self.par['MINSIZE']
      #deleted=list() #deleted id's

      for clumpId,pixels in clump.items():
         if len(pixels)<=minSize:
            #Set pixels as unusable
            for pos in pixels:
               caa[pos]=-1
            del clump[clumpId]
            #deleted.append(clumpId)
      if verbose: log.info('Number of remaining clumps: {0}'.format(len(clump)))

      ####refine 2, merging clumps
      """
      Amalgamate adjoining clumps if there is no significant dip between the
      clumps.
      """
      if verbose: log.info("[Stage 4] Merging clumps with low dip")

      """
      Get the minimum dip between two adjoining peaks necessary for the two
      peaks to be considered distinct.
      """
      minDip = self.par['MINDIP']*rms

      """
      Create a dictionary structure describing all the clumps boundaries in the
      supplied caa array. Then merge stage is applyed.
      """
      peaks,cols=self.clump_structs(clump,data,caa)
      clump,peaks,cols,caa=self.merge(clump,peaks,cols,caa,minDip)

      #ordering clumpId's (sequential id's)
      #verify if ordering is correct
      seqId=1
      for clumpId in clump:
         if clumpId!=seqId:
            clump[seqId]=clump.pop(clumpId)
            #update caa
            for pos in clump[seqId]:
               caa[pos]=seqId
         seqId+=1
      if verbose: log.info('Number of remaining clumps: {0}'.format(len(clump)))


      """
      Smooth the boundaries between the clumps. This cellular automata replaces
      each output pixels by the most commonly occuring value within a 3x3x3
      cube of input pixels centred on the output pixel. Repeat this process
      a number of times as given by configuration parameter CleanIter.
      """
      if verbose: log.info("[Stage 5] Smoothing boundaries")
      cleanIter = self.par['CLEANITER']
      for i in range(cleanIter):
         _caa = ca.smooth_boundary(caa)

         #positions where they were changes
         diff = caa!=_caa
         if not np.any(diff): del _caa; continue
         positions = map(tuple,np.array(np.where(diff)).T)

         #update clump struct
         for pos in positions:
            clump[caa[pos]].remove(pos)
            clump[_caa[pos]].append(pos)
         del caa
         caa = _caa
      return (caa,clump)



# Wrapper for a non-OO use of fellwalker
@support_nddata
def fellwalker(data,wcs=None,mask=None,unit=None,rms=0.0):
    fw=FellWalker()
    fw.defaultParams()
    clasar=fw.run(data,rms)
    return NDDataRef(clasar[0], uncertainty=None, mask=None, wcs=wcs, meta=None, unit=unit)
