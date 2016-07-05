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

