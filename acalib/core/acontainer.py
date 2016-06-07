from astropy import log

import acalib.io.formats as fmt

class AContainer:

    def __init__(self): 
        self.primary = None 
        self.adata = [] 
        self.atable = []

    def load(self,path): 
        fmt.load_to_cont(path,self)
    def save(self,path): 
        fmt.save_from_cont(path,self) 

