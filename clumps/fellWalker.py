import numpy as np
import copy
import matplotlib.pyplot as plt
import sys
from spectral import *


def get_params():
   maxJump=4
   minDip=3
   minSize=5
   return {'maxJump':maxJump,'minDip':minDip,'minSize':minSize}


def create_caa(data):
   caa=np.zeros_like(data)
   #Check invalid pixels (below threshold)
   shape=caa.shape
   rms=compute_rms(data)
   threshold=3*rms # here for now

   for i in range(shape[0]):
      for j in range(shape[1]):
         for k in range(shape[2]):
            if data[i,j,k]<threshold:
               caa[i,j,k]=-1
   print "CAA initialized successfully"
   return caa


def compute_rms(data):
   res=data[data < 0]
   fin=(res*res).sum()/len(res)
   return np.sqrt(fin)
       

def check_merge():
   return


def max_gradient(pos,data,caa):
   max_pos=pos
   max_val=data[pos]
   shape=data.shape

   for i in range(-1,2):
      for j in range(-1,2):
         for k in range(-1,2):
            neigh=(pos[0]+i,pos[1]+j,pos[2]+k) #position of neighbour point

            if i==j==k==0:
               #dont check it again
               continue
            elif neigh[0]<0 or neigh[1]<0 or neigh[2]<0 or neigh[0]>=shape[0] or neigh[1]>=shape[1] or neigh[2]>=shape[2]:
               #don't go outside of cube
               continue
            elif caa[neigh]==-1:
               #don't check unusable pixels
               continue
            elif data[neigh]>max_val:
               max_pos=neigh
               max_val=data[neigh]
   return max_pos

def verify_peak(pos,data,caa):
   max_pos=pos
   max_val=data[pos]
   shape=data.shape
   maxJump=get_params()['maxJump']

   for i in range(-maxJump,maxJump+1):
      for j in range(-maxJump,maxJump+1):
         for k in range(-maxJump,maxJump+1):
            neigh=(pos[0]+i,pos[1]+j,pos[2]+k) #position of neighbour point

            if abs(i)<=1 and abs(j)<=1 and abs(k)<=1:
               #don't check it again
               continue
            elif neigh[0]<0 or neigh[1]<0 or neigh[2]<0 or neigh[0]>=shape[0] or neigh[1]>=shape[1] or neigh[2]>=shape[2]:
               #don't go outside of cube
               continue
            elif caa[neigh]==-1:
               #don't check unusable pixels
               continue
            elif data[neigh]>max_val:
               max_pos=neigh
               max_val=data[neigh]
   return max_pos

def walkup(pos,path,data,caa):
   next_pos=max_gradient(pos,data,caa)

   if caa[next_pos]>=1:
      #Another ascent path reached
      path.append(next_pos)
      return path

   elif next_pos!=pos:
      #Keep walking up
      path.append(next_pos)
      return walkup(next_pos,path,data,caa)

   else:
      #Local peak reached
      local_max=next_pos
      new_max=verify_peak(local_max,data,caa)
      if local_max==new_max:
         #Significant peak reached
         return path
      else:
         #Just a noise peak, keep walking up   
         return walkup(new_max,path,data,caa)


def fellWalker(orig_cube):
   cube=copy.deepcopy(orig_cube)
   data=cube.data
   caa=create_caa(data)
   clump=dict()
   shape=data.shape
   top_id=0 #top clump id

   for i in range(shape[0]):
      for j in range(shape[1]):
         for k in range(shape[2]):
            # If unusable or already assigned, skip it
            if caa[i,j,k]!=0:
               continue
            path=list() # Ascent path pixels positions
            pos=(i,j,k)
            path.append(pos)
            path=walkup(pos,path,data,caa)
            print "nueva path:",path

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
            print "top_id",top_id

   #refine 1, removing small clumps
   minSize=get_params()['minSize']
   available_id=0 #trick to don't split ids

   for clumpId,pixels in clump.items():
      if len(pixels)<=minSize:
         #Set pixels as unusable
         for pos in pixels:
            caa[pos]=-1
         del clump[clumpId]

   return caa,clump