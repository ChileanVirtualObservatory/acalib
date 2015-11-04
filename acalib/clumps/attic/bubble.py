import numpy as np
from numpy.linalg import *
#from statistics import Gaussian
from scipy.cluster.vq import *
import copy
import matplotlib.pyplot as plt
import sys
from spectral import *
import matplotlib.cm as cm
#from sympy import symbols,cos,sin,exp,lambdify,diff
# Model: a,b,x0,y0,v0,dvx,dvy,sx,sy,sv,phi


def to_gauss((a,b,x0,y0,v0,phi,sx,sy,sv,dvx,dvy)):
   sphi2=np.square(np.sin(phi))
   cphi2=np.square(np.cos(phi))
   s2phi=np.sin(2*phi)
   sx2=np.square(sx)
   sy2=np.square(sy)
   sv2=np.square(sv)
   La=cphi2/sx2 + sphi2/sy2 + np.square(dvx)/sv2
   Ld=sphi2/sx2 + cphi2/sy2 + np.square(dvy)/sv2
   Lb=-s2phi/(2*sx2) + s2phi/(2*sy2) + dvx*dvy/sv2
   Lc=-dvx/sv2
   Le=-dvy/sv2
   Lf=1.0/sv2
   L=np.array([[La,Lb,Lc],[Lb,Ld,Le],[Lc,Le,Lf]])
   mu=[x0,y0,v0]
   return [a,b,mu,L]


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

def compute_rms(data):
   res=data[data < 0]
   fin=(res*res).sum()/len(res)
   return np.sqrt(fin)

def bubble_clump(orig_cube,size):
   plt.ion()
   plt.clf()
   cube=copy.deepcopy(orig_cube)
   std_pars=cube.standarize()
   syn=copy.deepcopy(cube)
   syn.data=np.empty_like(cube.data)
   rcube=copy.deepcopy(cube)
   rcube.data=np.empty_like(cube.data)
   D=size*np.array([cube.ra_delta,cube.dec_delta,cube.nu_delta])
   W=D*2.0
   # Compile Bubble
   center=np.zeros(3)
   center[0]=cube.ra_axis[int(len(cube.ra_axis)/2)]
   center[1]=cube.dec_axis[int(len(cube.dec_axis)/2)]
   center[2]=cube.nu_axis[int(len(cube.nu_axis)/2)]
   res=np.array([1.0,0,center[0],center[1],center[2],0,D[0],D[1],D[2],0,0])
   bubble=to_gauss(res)
   vect=np.empty((0,3))
   (features,index) = cube.feature_space(center,W) 
   val_fit=discrete_gauss(features,bubble)
   fit_cube=cube_data_unravel(val_fit,index)
   print fit_cube.max()
   print fit_cube.min()
   print fit_cube.shape
   (value_max,feature_max) = cube.max()
   mval=cube.data.mean()
   a=100*mval
   i=0
   while value_max > 2.06*mval:
      if i%100==0:
         print i
         print value_max, '>', mval
         plt.clf()
         plt.subplot(3, 2, 1)
         plt.imshow(cube.stack(),aspect='auto',origin='lower')
         plt.subplot(3, 2, 3)
         plt.imshow(cube.stack(axis=1),aspect='auto',origin='lower')
         plt.subplot(3, 2, 5)
         plt.imshow(cube.stack(axis=2),aspect='auto',origin='lower')
         plt.subplot(3, 2, 2)
         plt.imshow(syn.stack(),aspect='auto',origin='lower')
         plt.subplot(3, 2, 4)
         plt.imshow(syn.stack(axis=1),aspect='auto',origin='lower')
         plt.subplot(3, 2, 6)
         plt.imshow(syn.stack(axis=2),aspect='auto',origin='lower')
         plt.show()
         plt.pause(0.001)
         #mval=cube.data.mean()
      i=i+1
      #print feature_max, value_max
      index=cube.compute_window(feature_max,W)
      cube.add(-a*fit_cube,index)
      syn.add(a*fit_cube,index)
      vect=np.vstack((vect,feature_max))
      (value_max,feature_max) = cube.max()
   print "Bubbles =",i
   L=bubble[3]
   #print L
   Sig=inv(L)
   kmax=30
   error=np.empty(kmax)
   for k in range(kmax):
      print "k",k+1
      (codebook,dist)=kmeans(vect,k+1)
      (clus,ddis)=vq(vect,codebook)
      e1=0
      e2=0
      e3=0
      for i in range(k+1):
         v=vect[clus==i]
         #print "c"+str(i)+" size="+str(v.size/3)
         if (v.size==0):
            continue
         b=a*float(v.size/3)
         for mv in v:
             e1=e1+gauss_eval(v.T,(a*a,0,mv,L/2.0)).sum()
         ccv=np.cov(v.T)#/float(v.size/3)
         LL=inv(ccv+2*Sig)
         LL2=2*(ccv+Sig)
         e3=e3+b*b/np.sqrt(det(2*np.pi*LL2))
         e2=e2+a*gauss_eval(v.T,(b*a,0,codebook[i],LL)).sum()
      print e3, " - ",2*e2, " + ", e1
      error[k]=e1 - 2*e2 + e3
   k=np.argmin(error)
   print "Final K=",k+1
   (codebook,dist)=kmeans(vect,k+1)
   (clus,ddis)=vq(vect,codebook)
   (features,index) = cube.feature_space(center,np.array([1,1,1]))
   print index
   print features.shape
   for i in range(k+1):
         v=vect[clus==i]
         if (v.size==0):
            continue
         b=1.0/float(v.size/3)
         ccv=np.cov(v.T)
         LL=inv(ccv+Sig)
         W=2*LL
         val_fit=discrete_gauss(features,(b,0,codebook[i],LL))
         fit_cube=cube_data_unravel(val_fit,index)
         rcube.add(fit_cube,index)
   #rcube.unstandarize(std_pars)
   plt.clf() 
   zext=[cube.ra_axis[0],cube.ra_axis[-1],cube.dec_axis[0],cube.dec_axis[-1]]
   yext=[cube.ra_axis[0],cube.ra_axis[-1],cube.nu_axis[0],cube.nu_axis[-1]]
   xext=[cube.dec_axis[0],cube.dec_axis[-1],cube.nu_axis[0],cube.nu_axis[-1]]
   plt.subplot(4, 4, 1)
   plt.imshow(orig_cube.stack(),aspect='auto',origin='lower',extent=zext)
   plt.subplot(4, 4, 5)
   plt.imshow(orig_cube.stack(axis=1),aspect='auto',origin='lower',extent=yext)
   plt.subplot(4, 4, 9)
   plt.imshow(orig_cube.stack(axis=2),aspect='auto',origin='lower',extent=xext)
   plt.subplot(4, 4, 4)
   plt.imshow(rcube.stack(),aspect='auto',origin='lower',extent=zext)
   plt.subplot(4, 4, 8)
   plt.imshow(rcube.stack(axis=1),aspect='auto',origin='lower',extent=yext)
   plt.subplot(4, 4, 12)
   plt.imshow(rcube.stack(axis=2),aspect='auto',origin='lower',extent=xext)
   plt.subplot(4, 4, 2)
   plt.imshow(syn.stack(),aspect='auto',origin='lower',extent=zext)
   plt.subplot(4, 4, 6)
   plt.imshow(syn.stack(axis=1),aspect='auto',origin='lower',extent=yext)
   plt.subplot(4, 4, 10)
   plt.imshow(syn.stack(axis=2),aspect='auto',origin='lower',extent=xext)
   plt.subplot(4, 4, 13)
   plt.plot(error)
   plt.subplot(4, 4, 3)
   plt.xlim(zext[0],zext[1])
   plt.ylim(zext[2],zext[3])
   colors = iter(cm.rainbow(np.linspace(0, 1, k+1)))
   for i in range(k+1):
      v=vect[clus==i]
      plt.scatter(v[:,0],v[:,1], color=next(colors))
   plt.subplot(4, 4, 7)
   plt.xlim(yext[0],yext[1])
   plt.ylim(yext[2],yext[3])
   colors = iter(cm.rainbow(np.linspace(0, 1, k+1)))
   for i in range(k+1):
      v=vect[clus==i]
      plt.scatter(v[:,0],v[:,2], color=next(colors))
   plt.subplot(4, 4, 11)
   plt.xlim(xext[0],xext[1])
   plt.ylim(xext[2],xext[3])
   colors = iter(cm.rainbow(np.linspace(0, 1, k+1)))
   for i in range(k+1):
      v=vect[clus==i]
      plt.scatter(v[:,1],v[:,2], color=next(colors))
   plt.show()
   plt.pause(100)
