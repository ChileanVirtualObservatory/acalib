import astropy.units as u
from astropy import log
import numpy as np

FWHM_TO_SIGMA = 1. / (8 * np.log(2))**0.5

def _check(par,default):
    if not isinstance(par, u.Quantity):
        log.warning("Value '"+str(par)+"' without units, assuming '"+default.name+"'")
        par=par*default
    return par

def to_deg(par):
    # No quantity -> assumes deg
    par=_check(par,u.deg)
    return par.to(u.deg)

def to_hz(par):
    # No quantity -> assumes Hz
    par=_check(par,u.Hz)
    return par.to(u.Hz,equivalencies=u.spectral())

def to_rad(par):
    # No quantity -> assumes rad
    par=_check(par,u.rad)
    return par.to(u.rad)

def to_m_s(par):
    # No quantity -> assumes m/s
    par=_check(par,u.m/u.s)
    return par.to(u.m/u.s)

def vel_to_freq(vel,freq,equiv):
    freq=to_hz(freq)
    vel=to_m_s(vel)
    return vel.to(u.Hz,equivalencies=equiv(freq))
   
def fwhm_to_sigma(fwhm):
    return FWHM_TO_SIGMA*fwhm

def sigma_to_fwhm(fwhm):
    return fwhm/FWHM_TO_SIGMA

def to_hz_deg(grad,freq,equiv):
    grad=_check(grad,u.Hz/u.deg)
    vel=vel_to_freq(grad*u.deg,freq,equiv)
    return vel/u.deg

def ndslice(ndd, lower, upper):
    """ 
    N-Dimensional slicing.
    
    Arguments:
        ndd   -- an astropy.nddata.NDDataArray object.
        lower -- n-dimensional point as an n-tuple.
        upper -- n-dimensional point as an n-tuple.
    
    Returns:
        A sliced astropy.nddata.NDDataArray object.
        
    """
    lower = lower if lower is not None else np.zeros(ndd.ndim)
    upper = upper if upper is not None else ndd.shape
    return ndd[[slice(min(a,b), max(a,b)+1) for a,b in zip(lower, upper)]]

def adjust_index(relative, origin):
    """
    Adjusts an index relative to a subarray to an absolute
    index in the superarray.
    
    Arguments:
        origin   -- an n-dimensional index of a point as an n-tuple.
                    It should be the origin from which the relative
                    index was computed.
        relative -- an n-dimensional index of a point as an n-tuple.
                    The index to be adjusted.
    
    Returns:
        The relative index adjusted to the superarray as an n-tuple.
    """
    return tuple(np.array(origin) + np.array(relative))