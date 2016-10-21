

class Data(ndd.NDData):
    """
    A generic represenation of astronomical n-dimensional data array. Extends NDData.
    
    """

    def get_slab(self,lower=None,upper=None):
       pass

    def get_matching_slab(self,flux,lower,upper):
       pass

    def get_spectral_velocities(self,fqi=None,restfrq=None):
       pass

    def get_axes_ranges(self,lower=None,upper=None):
       pass

    def get_index_mesh(self,lower=None,upper=None):
       pass
    
    def get_index_features(self,lower=None,upper=None)
       pass

    def get_world_features(self,wcs=None,lower=None,upper=None):
       pass
 
    def get_fov_to_index(self,center,window):
       pass

    def get_rms(self,mask=None):
       pass

#HERE

#    def __init__(self,data,uncertainty,wcs,meta,unit):
#        ndd.NDData.__init__(self,data, uncertainty=None, mask=np.isnan(data), wcs=wcs, meta=meta, unit=unit)
#        self.data = ma.masked_array(data, mask=np.isnan(data))
 
#    @property
#    def data(self):
#        return self._data
        

#    @data.setter
#    def data(self, value):
#        self._data = value
    
#    def empty_like(self):
#    	dat=np.zeros_like(self.data)
#    	cb=Data(dat,self.wcs,self.meta,self.unit)
#    	return cb

#    def copy(self):
#        return copy.deepcopy(self)

#    def estimate_rms(self):
#       mm=self.data * self.data
#       rms=np.sqrt(mm.sum()*1.0/mm.count())
#       return rms


