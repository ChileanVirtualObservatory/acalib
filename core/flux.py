import astropy.units as u
import numpy as np

FWHM_TO_SIGMA = 1. / (8 * np.log(2))**0.5

def clump_to_gauss(pos,std,angle,freq,fwhm,gradient,equiv=u.doppler_radio):
   # Support for dimensionless parameters
   # Operational units: GHz and deg, except for angle (to radians for the trigonometric functions)
   if not isinstance(pos, u.Quantity):
      pos=pos*u.deg
   if not isinstance(std, u.Quantity):
      std=std*u.deg
   if not isinstance(freq, u.Quantity):
      freq=freq*u.GHz
   if not isinstance(fwhm, u.Quantity):
      fwhm=fwhm*u.km/u.s
   if not isinstance(gradient, u.Quantity):
      gradient=gradient*u.km/(u.s*u.deg)
   if not isinstance(angle, u.Quantity):
      angle=angle*u.deg
   # Conversions to operational units
   pos=(pos.to(u.deg)).value
   std=(std.to(u.deg)).value
   angle=(angle.to(u.rad)).value
   freq=(freq.to(u.GHz,equivalencies=u.spectral())).value
   sigma_freq=FWHM_TO_SIGMA * fwhm.to(freq.unit,equivalencies=equiv(freq))
   sigma_freq=(sigma_freq.to(u.GHz)).value
   grad_freq=(gradient*pos.unit).to(freq.unit,equivalencies=equiv(freq))*pos.unit
   grad_freq=(grad_freq.to(u.GHz/u.deg)).value
   # Construct the precision Matrix!
   sphi=np.sin(angle)
   cphi=np.cos(angle)
   R=np.array([[cphi,-sphi,-grad_freq[0]],[sphi,cphim-grad_freq[1]],[0,0,1]])
   D=np.diag([1./std[0],1./std[1],1./sigma_freq])
   RD=R.dot(D)
   P=RD.dot(RD.T)
   mu=np.array([pos[0],pos[1],freq])
   return (mu,P)

# Create a gassian flux on a cube with compact support 
def create_gauss_flux(wcs,mu,P,peak,cutoff):
   pass

