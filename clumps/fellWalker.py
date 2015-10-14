import copy
import sys
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

      """ Specific FellWalker parameters """
      #Longest jump to a higher neighbouring pixel
      self.par['MAXJUMP'] = 5
      #Minimum dip between distinct peaks
      self.par['MINDIP'] = 3
      #Minimum size (in pixels) of clumps
      self.par['MINSIZE'] = 20
      #Lowest gradient which marks the start of the walk
      self.par['FLATSLOPE'] = 0.001
      #Data values less than this are at "sea level"
      self.par['SEALEVEL'] = 0.

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
      caa = np.zeros_like(data.data).astype(np.int)
      #Check invalid pixels (below threshold)
      rms = self.par["RMS"]
      threshold = self.par['THRESH']*rms
      #Aditionally, NaN valued pixels are set as unusable -> filled(1)
      mask = np.array((data<threshold).filled(1))
      caa[mask] = -1
      print "CAA initialized successfully"
      return caa


   def compute_rms(self, data):
      mask = (data<0).filled(0)
      res = data[mask]
      fin = (res*res).sum()/len(res)
      return np.sqrt(fin)


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
      between them. Quit when an interation fails to merge any more clumps.
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
            print clumpId
            print neighId
            del tmpCols[neighId]
            tmpCols.update(cols[neighId])
            del tmpCols[clumpId]
            cols[neighId]=tmpCols
            for tmpId in cols:
               if clumpId in cols[tmpId]:
                  cols[tmpId][neighId]=cols[tmpId].pop(clumpId)

            #delete clumpId from peaks and update new peak of merged clumps
            tmpPeak=peaks.pop(clumpId)
            peaks[neighId]=np.max(tmpPeak,peaks[neighId])

            #Update caa
            for pos in clump[clumpId]:
               caa[pos]=neighId

            #Merge pixels in clump dictionary
            clump[neighId]+=(clump.pop(clumpId))
            print "Merged clumps {0} and {1}".format(clumpId,neighId)

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



   """
   *  Description:
   *     This function contracts (erodes) or expands (dilates) the pixels
   *     which are marked as "on" in the supplied input array using a cellular
   *     automata, and returns the result in an output array of the same shape
   *     and size as the input array.
   *
   *     Each output pixel value is created in turn as follows: If the
   *     corresponding input value has value equal to or greater than "magic",
   *     then the output pixel is set equal to "magic". Otherwise, a
   *     3x3x3 cube (or a 3x3 square for 2D arrays) is defined within the input
   *     array which is centred on the position of the current output pixel. The
   *     fraction of pixels within this input cube which are flagged as "on" is
   *     then counted. If this fraction is larger than "thresh" then the output
   *     pixel value is set to "on", otherwise it is set to "off". Optionally,
   *     an additional requirement for an output pixel to be set on is that
   *     the corresponding input pixel must be on.
   """

   def remove_isolate(self, inp, thresh, on, off, centre):
       out = np.empty_like(inp)

       shape = inp.shape

       for ox in range(shape[0]):
           for oy in range(shape[1]):
               for oz in range(shape[2]):
                   
                   if centre and inp[ox,oy,oz] == off:
                       """
                       If the corresponding input pixel is off, then the output must be also be
                       off if "centre" is true.
                       """
                       out[ox,oy,oz] = off

                   else:
                       """
                       Otherwise, loop round all input pixels in the neighbourhood of the current
                       output pixel, this is a cube of 3x3x3 input pixels, centred on the current
                       output pixel. Count how many of these input pixels are set to "on". If
                       the current output pixel is close to an edge of the array, there will be
                       fewer than 3x3x3 pixels in the cube. Count the total number of pixels
                       in the cube.
                       """
                       s = 0
                       tot = 0
                       for ix in range(ox-1, ox+2):
                           if ix < 0 or ix >= shape[0] : continue
                           for iy in range(oy-1, oy+2):
                               if iy < 0 or iy >= shape[1]: continue
                               for iz in range(oz-1, oz+2):
                                   if iz <0 or iz >= shape[2]: continue
                                   tot += 1
                                   if inp[ix,iy,iz] == on: s += 1
                       if np.float(s)/np.float(tot) > thresh:
                           out[ox,oy,oz] = on
                       else:
                           out[ox,oy,oz] = off
       return out


   def smooth_boundary(self, inp, clump):
       out = np.empty_like(inp)
       shape = inp.shape

       for ox in range(shape[0]):
           for oy in range(shape[1]):
               for oz in range(shape[2]):

                   party = dict()
                   maxvotes = 0
                   reached = False

                   for ix in range(ox-1, ox+2):
                       if ix < 0 or ix >= shape[0] : continue
                       for iy in range(oy-1, oy+2):
                           if iy < 0 or iy >= shape[1]: continue
                           for iz in range(oz-1, oz+2):
                               if iz < 0 or iz >= shape[2]: continue

                               inp_value = inp[ix, iy, iz]
                               #update party
                               if inp_value not in party:
                                   party[inp_value] = 1
                               else:
                                   party[inp_value] += 1

                               #update winner up to now
                               if party[inp_value] > maxvotes:
                                   winner = inp_value
                                   maxvotes = party[inp_value]

                               #verify if winner has been reached
                               if party[inp_value] > 13:
                                   winner = inp_value
                                   reached = True
                                   continue
                           if reached: continue
                       if reached: continue

                   #assign the winner
                   prev_clump = out[ox, oy, oz]
                   out[ox, oy, oz] = winner

                   #update clump dictionary
                   if prev_clump > 0:
                       clump[prev_clump].remove((ox,oy,oz))
                   if winner > 0:
                       clump[winner].append((ox,oy,oz))
       return out



   def fit(self, orig_cube):
      cube = orig_cube.copy()
      syn = orig_cube.empty_like()

      # Set the RMS, or automatically find an estimate for it
      if not self.par.has_key('RMS'):
         rms = cube.estimate_rms()
         self.par['RMS'] = rms

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
      Use a cellular automata to remove small isolated groups of usable
      pixels from the above "ipa" array. This allocates new memory for the
      cleaned up array.
      """
      print "Removing isolate clumps"
      frac = self.par['FRAC']
      on = self.par['ON']
      off = self.par['OFF']
      centre = self.par['CENTRE']
      caa = self.remove_isolate(caa, frac, on, off, centre)

      #Some constants
      rms = self.par['RMS']
      noise = self.par['THRESH']*rms
      flatSlope = self.par['FLATSLOPE']
      seaLevel = self.par['SEALEVEL']


      clump=dict()
      shape=data.shape
      top_id=0 #top clump id

      """
      Scan through the caa array, looking for usable pixels which have not
      yet been assigned to a clump (i.e. have a value of zero in caa).
      """
      print "Scaning pixel!"
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

               #Print some info
               print "id={0}, path: {1}".format(top_id,path)

      ###refine 1, removing small clumps
      print "\n Removing small clumps stage"
      minSize = self.par['MINSIZE']
      deleted=list() #deleted id's

      for clumpId,pixels in clump.items():
         if len(pixels)<=minSize:
            #Set pixels as unusable
            for pos in pixels:
               caa[pos]=-1
            del clump[clumpId]
            deleted.append(clumpId)
            print "removed clump {0}".format(clumpId)
         elif deleted:
            #If deleted have any index
            clump[deleted[0]]=clump.pop(clumpId)
            #update caa
            for pos in clump[deleted[0]]:
               caa[pos]=deleted[0]
            del deleted[0]
            deleted.append(clumpId)

      
      ####refine 2, merging clumps
      """
      Amalgamate adjoining clumps if there is no significant dip between the
      clumps.
      """
      print "\n Merge Stage"

      """
      Get the minimum dip between two adjoining peaks necessary for the two
      peaks to be considered distinct.
      """
      minDip = self.par['MINDIP']

      """
      Create a dictionary structure describing all the clumps boundaries in the
      supplied caa array. Then merge stage is applyed.
      """
      peaks,cols=self.clump_structs(clump,data,caa)
      clump,peaks,cols,caa=self.merge(clump,peaks,cols,caa,minDip)

      #ordering clumpId's (sequential id's)
      seqId=1
      for clumpId in clump:
         if clumpId!=seqId:
            clump[seqId]=clump.pop(clumpId)
            #update caa
            for pos in clump[seqId]:
               caa[pos]=seqId
         seqId+=1

      """
      Smooth the boundaries between the clumps. This cellular automata replaces
      each output pixels by the most commonly occuring value within a 3x3x3
      cube of input pixels centred on the output pixel. Repeat this process
      a number of times as given by configuration parameter CleanIter.
      """
      cleanIter = self.par['CLEANITER']
      for i in range(cleanIter): 
         caa = self.smooth_boundary(caa,clump)

      ####some statistics
      nclump=len(clump)
      print "\n SOME USEFUL DATA"
      print "Number of clumps:",nclump
      for clumpId,pixels in clump.items():
         print "Clump {0} has {1} pixels".format(clumpId,len(pixels))

      """
      Visualization of results
      """
      for pixels in clump.values():
         for pos in pixels:
            syn.data[pos]=data[pos]
            cube.data[pos]

      return caa,clump