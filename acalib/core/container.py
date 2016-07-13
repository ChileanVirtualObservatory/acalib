from astropy import log

from acalib.io.fits import save_fits_from_cont,load_fits_to_cont



class Container:
    """ Data structure that contains a list of NDData and astropy tables.
    
    It can load from and save to fits file format. To write a fits, please put in "primary"
    the main object, and in nndata and table all the extensions.
    """
    def __init__(self): 
        self.primary = None
        """Primary object""" 
        self.images = [] 
        """List of NDData object""" 
        self.tables = []
        """List of astropy tables""" 

    def load_fits(self,path): 
        load_fits_to_cont(path,self)
    def save_fits(self,path): 
        save_fits_from_cont(path,self) 


def load_fits(path):
    cont=Container()
    cont.load_fits(path)
    return cont
