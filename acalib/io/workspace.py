
from astropy.io import fits
from astropy.nddata import NDData
from astropy.table import Table
from astropy.wcs import wcs
import numpy as np
from astropy import log
from astropy import __version__
import time
import tempfile
import formats 

def create(name):
   ws=dict()
   ws['workspace']=name
   return ws

_ws_df=create("DEFAULT")

def elements(ws=_ws_df):
   retval=ws.copy()
   del retval['workspace']
   return retval

def import_file(path,ws=_ws_df):
   formats.load_to_ws(path,ws)

def real_dims(ndd):
   shape=[]
   dim=0

   for i in range(ndd.data.ndim):
      if ndd.data.shape[i] != 1:
         shape.append(ndd.data.shape[i])
         dim+=1
   if dim==1:
      otype="Spectra"
   elif dim==2:
      otype="Image"
   elif dim>=3:
      otype="Cube"
   else:
      log.warning("NDData of 0 dimension? ignoring...")
      return
   return (dim,shape,otype)

