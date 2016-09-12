import astropy.units as u
from astropy import log
import numpy as np

FWHM_TO_SIGMA = 1. / (8 * np.log(2))**0.5

#TODO Document

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

def gclump_to_wcsgauss(pos,std,angle,freq,fwhm,gradient,equiv=u.doppler_radio):
   # Parameter sanitization
   pos=to_deg(pos)
   std=to_deg(std)
   angle=to_rad(angle)
   freq=to_hz(freq)
   #print "fwhm",fwhm
   #print "freq",freq
   sigma=fwhm_to_sigma(freq - vel_to_freq(fwhm,freq,equiv))
   #print "sigma",sigma
   grad= freq/u.deg -  to_hz_deg(gradient,freq,equiv)
   # get Values
   pos=pos.value
   std=std.value
   angle=angle.value
   freq=freq.value
   sigma=sigma.value
   grad=grad.value
   # Construct the precision Matrix!
   sphi=np.sin(angle)
   cphi=np.cos(angle)
   R=np.array([[cphi,-sphi,-grad[0]],[sphi,cphi,-grad[1]],[0,0,1]])
   D=np.diag([1./std[0],1./std[1],1./sigma])
   RD=R.dot(D)
   P=RD.dot(RD.T)
   mu=np.array([pos[0],pos[1],freq])
   return (mu,P)


