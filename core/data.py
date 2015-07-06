import numpy as np
import copy

from astropy import constants as const
from astropy import units as u
from astropy.io import fits 
from astropy import log
import astropy.nddata as ndd
import numpy.ma as ma
import astropy.wcs as astrowcs
import matplotlib.pyplot as plt
import parameter as par
import time
import scipy.ndimage

class AcaData(ndd.NDData):
    """
    A generic represenation of astronomical data.
    A spectra is a 3D cube with ra_axis and dec_axis of size 1
    An image is a 3D cube with nu_axis of size 1
    A spectroscopic cube is a 3D cube
    Stokes 4D cubes are not supported.
    """
    def __init__(self,data,wcs,meta,unit):
        ndd.NDData.__init__(self,data, uncertainty=None, mask=np.isnan(data), wcs=wcs, meta=meta, unit=unit)
        self.data = ma.masked_array(data, mask=np.isnan(data))
 
    @property
    def data(self):
        return self._data
        

    @data.setter
    def data(self, value):
        self._data = value
    
#    def _init_from_params(self, pos, ang_res, fov, freq, spe_res, bw):
#        # Create a new WCS object.
#        pos=par.to_deg(pos)
#        ang_res=par.to_deg(ang_res)
#        fov=par.to_deg(fov)
#        #print pos,ang_res,fov
#        freq=par.to_hz(freq)
#        spe_res=par.to_hz(spe_res)
#        bw=par.to_hz(bw)
#        #print freq,spe_res,bw
#        wcs = astrowcs.WCS(naxis=3)
#        wcs.wcs.crval = np.array([pos[0].value, pos[1].value,freq.value])
#        wcs.wcs.cdelt = np.array([ang_res[0].value, ang_res[1].value,spe_res.value])
#        mm = np.array([int(abs(fov[0]/ang_res[0])),int(abs(fov[1]/ang_res[1])),int(abs(bw/spe_res))])
#        wcs.wcs.crpix = mm/2.0
#        wcs.wcs.ctype = ["RA---SIN", "DEC--SIN","FREQ"]
#        data=np.zeros((mm[2],mm[1],mm[0]))
#        #wcs.wcs.print_contents()
#        ndd.NDData.__init__(self,data, uncertainty=None, mask=np.isnan(data), wcs=wcs, meta=None, unit=u.Jy/u.beam)

    def get_wcs_limits(self,axis):
       lower=self.wcs.wcs_pix2world([[0,0,0]], 0) - self.wcs.wcs.cdelt/2.0
       shape=self.data.shape
       shape=[shape[::-1]]
       upper=self.wcs.wcs_pix2world(shape, 1) + self.wcs.wcs.cdelt/2.0
       #print lower,upper
       return (lower[0][axis],upper[0][axis])    

    def copy(self):
        return copy.deepcopy(self)

    # TODO: get_flux should be for any slice, need more parameters
    def get_flux(self):
    	return np.sum(self.data)   
    
    
    def empty_like(self):
    	dat=np.zeros_like(self.data)
    	cb=AcaData(dat,self.meta)
    	return cb
    
    def scale(self, scale):
    	zflag = 0
    	start_time = time.time()
    	if (scale == 1):
    		return self.data
    	elif (scale < 1):
    		new_data = np.zeros((round(len(self.data)*scale+1),round(len(self.data[0])*scale+1), round(len(self.data[0][0])*scale+1)))          
    		for z in self.data:
    			yflag = 0
    			for y in z:
    				xflag = 0
    				for x in y:
    					new_data[round(zflag*scale)][round(yflag*scale)][round(xflag*scale)] = x
    					xflag+=1
    				yflag+=1
    			zflag+=1				
    		#for z in range(len(self.data)):
    		#	for y in range(len(self.data[0])):
    		#			for x in range(len(self.data[0][0])):
    		#					new_data[round(z*scale)][round(y*scale)][round(x*scale)] = self.data[z][y][x]
    		print("--- %s seconds ---" % (time.time() - start_time))
    		return (new_data/np.sum(new_data))*np.sum(self.data)
    	else:
    		new_data = np.zeros((round(len(self.data)*scale),round(len(self.data[0])*scale), round(len(self.data[0][0])*scale)))          
    		for z in self.data:
    			new_data[zflag] = scipy.ndimage.zoom(z,round(scale),order=3)
    			zflag+=1
    		print("--- %s seconds ---" % (time.time() - start_time))
    		return new_data     
    							 
    	 
    def _slice(self,lower,upper):
    		if lower==None:
    				lower=(0,0,0)
    		if upper==None:
    				upper=self.data.shape
    		if isinstance(lower,tuple):
    				lower=np.array(lower)
    		if isinstance(upper,tuple):
    				upper=np.array(upper)
    		llc=lower < 0
    		ulc=lower > self.data.shape
    		luc=upper < 0
    		uuc=upper > self.data.shape
    		if llc.any():
    				log.warning("Negative lower index "+str(lower)+". Correcting to zero.")
    				lower[llc]=0
    		if ulc.any():
    				log.warning("Lower index out of bounds "+str(lower)+" > "+str(self.data.shape)+". Correcting to max.")
    				upper[ulc]=self.data.shape[ulc]
    		if luc.any():
    				log.warning("Negative upper index "+str(upper)+". Correcting to zero.")
    				lower[luc]=0
    		if uuc.any():
    				log.warning("Upper index out of bounds "+str(upper)+" > "+str(self.data.shape)+". Correcting to max.")
    				upper[uuc]=np.array(self.data.shape)[uuc]
    		return [slice(lower[0],upper[0]),slice(lower[1],upper[1]),slice(lower[2],upper[2])]
    			
    def get_stacked(self,lower=None,upper=None,axis=(0)):
                sli=self._slice(lower,upper)
                # TODO: nan values must be excluded using mask (data.sum), but they are not working!
    		return np.nansum(self.data[sli],axis=axis)
    
    def add_flux(self,flux,lower=None,upper=None):
    		sli=self._slice(lower,upper)
    		fl=np.array([0,0,0])
    		fu=np.array(flux.shape)
    		for i in range(0,3):
    			 if sli[i].start == 0:
    					fl[i]=flux.shape[i] - sli[i].stop
    			 if sli[i].stop == self.data.shape[i]:
    					fu[i]=sli[i].stop - sli[i].start
    		self.data[sli]+=flux[fl[0]:fu[0],fl[1]:fu[1],fl[2]:fu[2]]
    
    def max(self):
    		# TODO: here we should use only self.data.argmax(), but nanargmax is used
    		# because the bloody masked arrays are not working in the NDData
    		index=np.unravel_index(np.nanargmax(self.data),self.data.shape)
    		y=self.data[index]
    		return (y,index)
    
    def min(self):
    		# TODO: here we should use only self.data.argmin(), but nanargmin is used
    		# because the bloody masked arrays are not working in the NDData
    		index=np.unravel_index(np.nanargmin(self.data),self.data.shape)
    		y=self.data[index]
    		return (y,index)
    
    def index_to_wcs(self,index):
    		val=self.wcs.wcs_pix2world([index[::-1]],0)
    		if val.shape[0]==1: val=val[0]
    		return val
    
    def get_axis_names(self):
    		return self.wcs.axis_type_names
    
    def get_features(self,lower=None,upper=None):
    		sli=self._slice(lower,upper)
    		x=np.arange(sli[0].start,sli[0].stop)
    		y=np.arange(sli[1].start,sli[1].stop)
    		z=np.arange(sli[2].start,sli[2].stop)
    		xyz=np.meshgrid(x,y,z,indexing='ij')
    		ii=np.empty((3,len(x)*len(y)*len(z)))
    		ii[2]=xyz[0].ravel()
    		ii[1]=xyz[1].ravel()
    		ii[0]=xyz[2].ravel()
    		f=self.wcs.wcs_pix2world(ii.T,0)
    		f=f.T
    		return f
    
    def get_slice(self,lower=None,upper=None):
    		sli=self._slice(lower,upper)
    		return np.array([self.data[sli[0],sli[1],sli[2]]])
    
    def index_from_window(self,wcs_center,wcs_window):
		   ld=np.rint(self.wcs.wcs_world2pix([wcs_center-wcs_window],0))
		   lu=np.rint(self.wcs.wcs_world2pix([wcs_center+wcs_window],0))
		   lower=np.array([ld,lu]).min(axis=0)
		   upper=np.array([ld,lu]).max(axis=0)
		   return (lower[0][::-1],upper[0][::-1])
    
#    def _add_HDU(self, hdu):
#    		self.hdulist.append(hdu)
    
#    def save_fits(self, filename):
#    		""" Simple as that... saves the whole cube """
#    		# TODO: Check this, I think we should add STOKES to be 100% compatible to ALMA
#    		# TODO: Add a proper wcs._to_header or wcs._to_fits...
#    		self.hdulist.writeto(filename, clobber=True)
#    
#    def _updatefig(self, j):
#    		""" Animate helper function """
#    		self.im.set_array(self.data[j, :, :])
#    		return self.im,
    
#    def animate(self, inte, rep=True):
#    		#TODO: this is not ported to the new wcs usage: maybe we must use wcsaxes to plot the wcs information...
#    		""" Simple animation of the cube.
#    				- inte       : time interval between frames
#    				- rep[=True] : boolean to repeat the animation
#    			"""
#    		fig = plt.figure()
#    		self.im = plt.imshow(self.data[0, :, :], cmap=plt.get_cmap('jet'), vmin=self.data.min(), vmax=self.data.max(), \
#    												 extent=(
#    														 self.alpha_border[0], self.alpha_border[1], self.delta_border[0],
#    														 self.delta_border[1]))
#    		ani = animation.FuncAnimation(fig, self._updatefig, frames=range(len(self.freq_axis)), interval=inte, blit=True,
#    																	repeat=rep)
#    		plt.show()

#    def max_energy(self,sc,idx):
#          target=self.data[idx[4]:idx[5],idx[2]:idx[3],idx[0]:idx[1]]
#          if target.shape != sc.shape:
#            si=np.array([0,sc.shape[2],0,sc.shape[1],0,sc.shape[0]])
#            mm=target.min()
#            datum=mm*np.ones_like(sc)
#            if idx[4] == 0:
#               si[4]=si[5]  - idx[5]
#            if idx[2] == 0:
#               si[2]=si[3]  - idx[3]
#            if idx[0] == 0:
#               si[0]=si[1]  - idx[1]
#            if idx[5] == self.nu_axis.size:
#               si[5] =idx[5] - idx[4] 
#            if idx[3] == self.dec_axis.size:
#               si[3] =idx[3] - idx[2] 
#            if idx[1] == self.ra_axis.size:
#               si[1] =idx[1] - idx[0]
#            datum[si[4]:si[5],si[2]:si[3],si[0]:si[1]]=target
#            #max_energy=(self.data[idx[4]:idx[5],idx[2]:idx[3],idx[0]:idx[1]]/sc[si[4]:si[5],si[2]:si[3],si[0]:si[1]]).min()
#          else:
#            datum=target
#          max_energy=(datum/sc).min()
#          return max_energy

