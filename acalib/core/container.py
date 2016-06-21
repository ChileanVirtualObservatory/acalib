from astropy import log

from acalib.io.fits import save_fits_from_cont,load_fits_to_cont


class Container:

    def __init__(self): 
        self.primary = None 
        self.nddata = [] 
        self.table = []

    def load_fits(self,path): 
        load_fits_to_cont(path,self)
    def save_fits(self,path): 
        save_fits_from_cont(path,self) 

