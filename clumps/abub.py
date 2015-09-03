import numpy as np
from numpy.linalg import *
from scipy.cluster.vq import *
import matplotlib.pyplot as plt
from spectral import *
import matplotlib.cm as cm
import scipy.stats
import scipy.spatial.distance as dist
import scipy.cluster.hierarchy as hier

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

# Simple hill-climbing to find the nearby region with more energy
def improve_energy(a,index,bubble,cube,weight):
  mods=np.empty((0,6))
  sizes=np.array([index[1]-index[0],index[3]-index[2],index[5]-index[4]])
  if index[0] > 0:
     if index[1] == cube.ra_axis.size and sizes[0] < bubble.shape[2]:
        mods=np.vstack((mods,[-1,0,0,0,0,0]))
     else:
        mods=np.vstack((mods,[-1,-1,0,0,0,0]))
  if index[1] < cube.ra_axis.size:
     if index[0] == 0 and sizes[0] < bubble.shape[2]:
        mods=np.vstack((mods,[0,1,0,0,0,0]))
     else:
        mods=np.vstack((mods,[1,1,0,0,0,0]))
  if index[2] > 0 :
     if index[3] == cube.dec_axis.size and sizes[1] < bubble.shape[1]:
        mods=np.vstack((mods,[0,0,-1,0,0,0]))
     else:
        mods=np.vstack((mods,[0,0,-1,-1,0,0]))
  if index[3] < cube.dec_axis.size:
     if index[2] == 0 and sizes[1] < bubble.shape[1]:
        mods=np.vstack((mods,[0,0,0,1,0,0]))
     else:
        mods=np.vstack((mods,[0,0,1,1,0,0]))
  if index[4] > 0 :
     if index[5] == cube.nu_axis.size and sizes[2] < bubble.shape[0]:
        mods=np.vstack((mods,[0,0,0,0,-1,0]))
     else:
        mods=np.vstack((mods,[0,0,0,0,-1,-1]))
  if index[5] < cube.nu_axis.size:
     if index[4] == 0 and sizes[2] < bubble.shape[0]:
        mods=np.vstack((mods,[0,0,0,0,0,1]))
     else:
        mods=np.vstack((mods,[0,0,0,0,1,1]))
  newidx=np.array([])
  for m in mods:
     idx=index+m
     ap=cube.max_energy(bubble,idx)*weight
     if ap > a:
        a=ap
        newidx=idx
  if newidx.size != 0 :
     (a,index)=improve_energy(a,newidx,bubble,cube,weight)
  return (a,index)

# Find where to extract the next bubble
def next_bubble(cube,window,bubble,weight):
  (value_max,max_pos) = cube.max()
  index=np.array(cube.compute_window(max_pos,window))
  a=cube.max_energy(bubble,index)*weight
  (a,index)=improve_energy(a,index,bubble,cube,weight)
  return (a,index)

def plot_init():
   plt.ion()
   plt.clf()

def plot_iter_status(orig,cube,syn,vect,ener,varia,snr):
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
   plt.ylabel("cumsum")
   cs=ener.cumsum()
   plt.plot(cs)
   #plt.plot(ssnr*np.ones(ener.size),'r')
   plt.subplot(4, 4, 15)
   plt.xlabel("iter")
   plt.ylabel("variance")
   plt.plot(varia)
   plt.subplot(4, 4, 16)
   plt.xlabel("iter")
   plt.ylabel("SNR")
   plt.plot(snr)
   
   plt.show()
   plt.pause(0.001)


def bubble_fit(orig,size,weight=0.5,verbose=True,plot=True,report_every=100,plot_every=100):
   # Standarize 0-1 and make an empty cube:
   cube=orig.copy()
   std_pars=cube.standarize()
   syn=0
   snr=0
   if plot:
      plot_init()
      syn=cube.empty_like()
   
   if verbose:  
      print "Standarization Report:"
      print "--> Constants", std_pars

   # Create a generic 0-1 bubble:
   #   a. Select the middle pixel of the cube
   mu=np.zeros(3)
   mu[0]=cube.ra_axis[int(len(cube.ra_axis)/2)]
   mu[1]=cube.dec_axis[int(len(cube.dec_axis)/2)]
   mu[2]=cube.nu_axis[int(len(cube.nu_axis)/2)]
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
   varia=np.vstack((varia,cube.data.std()))
   snr=np.empty((0,1))
   
   i=0
   tot_a=0
   if verbose:
      print "Bubble Extraction Iteration:"
   while True:
      # Obtain the next bubble
      (a,index)=next_bubble(cube,window,bubble,weight)
      cube.add(-a*bubble,index)
      if plot:
         syn.add(a*bubble,index)
      vect=np.vstack((vect,cube.index_center(index)))
      ener=np.vstack((ener,a))
      tot_a+=a
      varia=np.vstack((varia,cube.data.std()))
      snr=np.vstack((snr,cube.snr()))
      
      # Text Reports
      if verbose and i%report_every==0:
         # Compute statistics
         totener=ener.sum()
         print "--> bubbles =",i
         print "    * total energy =",tot_a
         print "    * next energy =",a
         print "    * variance =",varia[i]
         print "    * snr =",snr[i]
         #print "    * variance =",varia[i/report_every]
         #print "    * snr =",snr[i/report_every]
      if plot and i%plot_every==0:
         plot_iter_status(orig,cube,syn,vect,ener,varia,snr)
      
      # End Iteration
      i=i+1
   # Heriarchical Clustering
   # Compute the condensated eucledian distance matrix
   M=dist.pdist(vect)
   # Compute the weighting factors of the distance matrix
   W=dist.pdist(ener,lambda u,v: u+v)
   # Perform agglomerative clustering
   Z=hier.linkage(M*W)
   #hier.dendrogram(Z)
   zext=[cube.ra_axis[0],cube.ra_axis[-1],cube.dec_axis[0],cube.dec_axis[-1]]
   yext=[cube.ra_axis[0],cube.ra_axis[-1],cube.nu_axis[0],cube.nu_axis[-1]]
   xext=[cube.dec_axis[0],cube.dec_axis[-1],cube.nu_axis[0],cube.nu_axis[-1]]
   n=1
   ll=5
   offset=2
   if plot:
      plt.pause(10)
   plt.clf()
   for j in range(ll):
      k=j+offset
      T=hier.fcluster(Z,k,criterion='maxclust')
      plt.subplot(ll, 3, n)
      plt.xlim(zext[0],zext[1])
      plt.ylim(zext[2],zext[3])
      colors = iter(cm.rainbow(np.linspace(0, 1, k+1)))
      for i in range(k+1):
         v=vect[T==i]
         print v
         plt.scatter(v[:,0],v[:,1], color=next(colors))
      n=n+1
      plt.subplot(ll, 3, n)
      plt.xlim(yext[0],yext[1])
      plt.ylim(yext[2],yext[3])
      colors = iter(cm.rainbow(np.linspace(0, 1, k+1)))
      for i in range(k+1):
         v=vect[T==i]
         plt.scatter(v[:,0],v[:,2], color=next(colors))
      n=n+1
      plt.subplot(ll, 3, n)
      plt.xlim(xext[0],xext[1])
      plt.ylim(xext[2],xext[3])
      colors = iter(cm.rainbow(np.linspace(0, 1, k+1)))
      for i in range(k+1):
         v=vect[T==i]
         plt.scatter(v[:,1],v[:,2], color=next(colors))
      n=n+1
   plt.show()
   plt.pause(500)
   return (vect,ener)

