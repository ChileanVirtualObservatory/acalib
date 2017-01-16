from astropy.table import Table, Column
import numpy as np
from .axes import axes_names

def _pix_table_creator(values,wcs):
    tab = Table()
    names = axes_names(None, wcs)
    for i in range(names.size):
       tab[names[i]]=Column(values[:,i],unit=u.pix)
    return tab

def _world_table_creator(values,wcs):
    uvec=np.array(wcs.wcs.cunit)[::-1]
    values=np.fliplr(values)
    tab = Table()
    names = axes_names(None, wcs)
    for i in range(names.size):
       tab[names[i]]=Column(values[:,i],unit=uvec[i])
    return tab


def _unitize(vec,wcs):
    uvec=np.array(wcs.wcs.cunit)[::-1]
    return vec*uvec
