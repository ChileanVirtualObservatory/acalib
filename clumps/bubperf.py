import numpy as np
from numpy.linalg import *
from scipy.cluster.vq import *
import matplotlib.pyplot as plt
from spectral import *
import matplotlib.cm as cm
import scipy.stats

from joblib import Parallel, delayed  
import multiprocessing

def gauss_eval(features,(a,b,mu,L)):
   C=np.empty_like(features)
   C[0]=features[0] - mu[0]
   C[1]=features[1] - mu[1]
   C[2]=features[2] - mu[2]
   V=C*(L.dot(C))
   quad=V.sum(axis=0)
   v=np.exp(-quad/2.0)
   retval=b + a*v*np.sqrt(det(L)/np.power(2*np.pi,3));
   return retval


#TODO: Replace this for erf computation... or beam, or whatever...
def discrete_gauss(features,(a,b,mu,L)):
   C=np.empty_like(features)
   C[0]=features[0] - mu[0]
   C[1]=features[1] - mu[1]
   C[2]=features[2] - mu[2]
   V=C*(L.dot(C))
   quad=V.sum(axis=0)
   v=np.exp(-quad/2.0)
   v=v/v.sum()
   retval=b + a*v;
   return retval

def compute_energy(bubble,cube,pos,energies):
  index=index_by_pos(pos,bubble,cube)
  energies[pos[0],pos[1],pos[2]]= cube.max_energy(bubble,index)

def update_energies(energies,cube,bubble,pos,region=np.array([]),verbose=False):
  num_cores = multiprocessing.cpu_count()
  if region.size == 0:
     # Compute region
     region=[0]*6
     region[0]=int(max(0,pos[2]-bubble.shape[2]))
     region[1]=int(min(cube.ra_axis.size,pos[2]+bubble.shape[2]))
     region[2]=int(max(0,pos[1]-bubble.shape[1]))
     region[3]=int(min(cube.dec_axis.size,pos[1]+bubble.shape[1]))
     region[4]=int(max(0,pos[0]-bubble.shape[0]))
     region[5]=int(min(cube.nu_axis.size,pos[0]+bubble.shape[0]))
  XX,YY,ZZ=np.meshgrid(range(region[0],region[1]),range(region[2],region[3]),range(region[4],region[5]))
  PP=np.array([ZZ.flatten(),YY.flatten(),XX.flatten()]).T
  i=0
  ss=PP.shape[0]
  for p in PP:
     i=i+1
     if i%100000==0:
       print i*100/ss,"%"
     compute_energy(bubble,cube,p,energies)
#for x in range(region[0],region[1]):
#   if verbose:
#      print (x-region[0])*100/(region[1]-region[0]),"%"
#   for y in range(region[2],region[3]):
#      
#      for z in range(region[4],region[5]):
#         pos=np.array([z,y,x])
#if x==0:
#   print index
#if not verbose:
#   print index,cube.max_energy(bubble,index)
#  Parallel(n_jobs=num_cores)(delayed(compute_energy)(bubble,cube,pos,energies) for pos in PP)

def index_by_pos(pos,bubble,cube):
  index=np.empty(6)
  index[0]=np.maximum(0,pos[2]-bubble.shape[2]/2)
  index[1]=np.minimum(cube.ra_axis.size,pos[2]+bubble.shape[2]/2+1)
  index[2]=np.maximum(0,pos[1]-bubble.shape[1]/2)
  index[3]=np.minimum(cube.dec_axis.size,pos[1]+bubble.shape[1]/2+1)
  index[4]=np.maximum(0,pos[0]-bubble.shape[0]/2)
  index[5]=np.minimum(cube.nu_axis.size,pos[0]+bubble.shape[0]/2+1)
  return index

# Find where to extract the next bubble
def next_bubble(cube,bubble,weight,energies):
  pos=np.unravel_index(energies.argmax(),energies.shape)
  index=index_by_pos(pos,bubble,cube)
  a=cube.max_energy(bubble,index)*weight
  return (a,index,pos)

def plot_init():
   plt.ion()
   plt.clf()

def plot_iter_status(orig,cube,syn,vect,ener,varia,entro):
   plt.clf()
   # First
   plt.subplot(4, 4, 1)
   plt.xlabel("ra")
   plt.ylabel("dec")
   plt.imshow(orig.stack(),aspect='auto',origin='lower')
   plt.subplot(4, 4, 2)
   plt.xlabel("ra")
   plt.ylabel("dec")
   plt.imshow(cube.stack(),aspect='auto',origin='lower')
   plt.subplot(4, 4, 3)
   plt.xlabel("ra")
   plt.ylabel("dec")
   plt.imshow(syn.stack(),aspect='auto',origin='lower')
   plt.subplot(4, 4, 4)
   plt.xlabel("ra")
   plt.ylabel("dec")
   plt.xlim(0,1)
   plt.ylim(0,1)
   plt.scatter(vect[:,0],vect[:,1])
   # Second
   plt.subplot(4, 4, 5)
   plt.xlabel("ra")
   plt.ylabel("nu")
   plt.imshow(orig.stack(axis=1),aspect='auto',origin='lower')
   plt.subplot(4, 4, 6)
   plt.xlabel("ra")
   plt.ylabel("nu")
   plt.imshow(cube.stack(axis=1),aspect='auto',origin='lower')
   plt.subplot(4, 4, 7)
   plt.xlabel("ra")
   plt.ylabel("nu")
   plt.imshow(syn.stack(axis=1),aspect='auto',origin='lower')
   plt.subplot(4, 4, 8)
   plt.xlabel("ra")
   plt.ylabel("nu")
   plt.xlim(0,1)
   plt.ylim(0,1)
   plt.scatter(vect[:,0],vect[:,2])
   # Third
   plt.subplot(4, 4, 9)
   plt.xlabel("nu")
   plt.ylabel("dec")
   plt.imshow(orig.stack(axis=2).T,aspect='auto',origin='lower')
   plt.subplot(4, 4, 10)
   plt.xlabel("nu")
   plt.ylabel("dec")
   plt.imshow(cube.stack(axis=2).T,aspect='auto',origin='lower')
   plt.subplot(4, 4, 11)
   plt.xlabel("nu")
   plt.ylabel("dec")
   plt.imshow(syn.stack(axis=2).T,aspect='auto',origin='lower')
   plt.subplot(4, 4, 12)
   plt.xlabel("nu")
   plt.ylabel("dec")
   plt.xlim(0,1)
   plt.ylim(0,1)
   plt.scatter(vect[:,2],vect[:,1])
   plt.subplot(4, 4, 13)
   plt.xlabel("iter")
   plt.ylabel("energy")
   plt.plot(ener)
   plt.subplot(4, 4, 14)
   plt.xlabel("iter")
   plt.ylabel("cum(energy)")
   plt.plot(ener.cumsum())
   plt.subplot(4, 4, 15)
   plt.xlabel("iter")
   plt.ylabel("variance")
   plt.plot(varia)
   plt.subplot(4, 4, 16)
   plt.xlabel("iter")
   plt.ylabel("entropy")
   plt.plot(entro)
   
   plt.show()
   plt.pause(0.001)


def bubble_perfect(orig,size,weight=0.5,verbose=True,plot=True,report_every=100,plot_every=100):
   

   
   if plot:
      plot_init()
   
   # Standarize 0-1 and make an empty cube:
   cube=orig.copy()
   std_pars=cube.standarize()
   syn=cube.empty_like()
   eee=cube.empty_like()
   
   if verbose:  
      print "Standarization Report:"
      print "--> Constants", std_pars

   # Create a generic 0-1 bubble:
   #   a. Select the middle pixel of the cube
   pos=np.zeros(3)
   mu=np.zeros(3)
   pos[2]=int(len(cube.ra_axis)/2)
   pos[1]=int(len(cube.dec_axis)/2)
   pos[0]=int(len(cube.nu_axis)/2) 
   mu[0]=cube.ra_axis[pos[2]]
   mu[1]=cube.dec_axis[pos[1]]
   mu[2]= cube.nu_axis[pos[0]]
   #   b. define the bubble size and window size
   resol  = np.array([cube.ra_delta,cube.dec_delta,cube.nu_delta])
   sigmas = size*resol
   window = 2.0*sigmas
   #   c. create the inverse of the covariance matrix (precision matrix)
   lambd  = np.array([[1.0/np.square(sigmas[0]),0,0],[0,1.0/np.square(sigmas[1]),0],[0,0,1.0/np.square(sigmas[2])]])
   #   d. select the feature space and proyect the gaussian bubble there
   (features,index) = cube.feature_space(mu,window) 
   feat_bubble=discrete_gauss(features,(1,0,mu,lambd))
   #   e. unravel the features to obtain the bubble in a cube
   bubble=cube_data_unravel(feat_bubble,index)
   
   if verbose:  
      print "Bubble Construction Report:"
      print "--> Sigmas", sigmas
      print "--> Pixels", bubble.shape

   # Create empty vectors for stacking positions and energies
   vect=np.empty((0,3))
   ener=np.empty((0,1))

   # Create empty vectors for stacking statistics
   varia=np.empty((0,1))
   entro=np.empty((0,1))
   energies=np.empty_like(cube.data)
   region=np.array([0,cube.ra_axis.size,0,cube.dec_axis.size,0,cube.nu_axis.size])
   update_energies(energies,cube,bubble,pos,region,verbose)
   i=0
   if verbose:
      print "Bubble Extraction Iteration:"
   while True :
      # Obtain the next bubble
      (a,index,pos)=next_bubble(cube,bubble,weight,energies)
      cube.add(-a*bubble,index)
      syn.add(a*bubble,index)
      update_energies(energies,cube,bubble,pos)
      vect=np.vstack((vect,cube.index_center(index)))
      ener=np.vstack((ener,a))

      
      # Text Reports
      if verbose and i%report_every==0:
         # Compute statistics
         totener=ener.sum()
         varia=np.vstack((varia,(cube.data/(1.0-totener)).std()))
         entro=np.vstack((entro,scipy.stats.entropy( (cube.data/(1.0-totener)).flatten())))
         print "--> bubbles =",i
         print "    * next energy =",a
         print "    * variance =",varia[i/report_every]
         print "    * entropy =",entro[i/report_every]
      if plot and i%plot_every==0:
         plot_iter_status(orig,cube,syn,vect,ener,varia,entro)
      
      # End Iteration
      i=i+1


