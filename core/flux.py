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
   R=np.array([[cphi,-sphi,-grad_freq[0]],[sphi,cphi,-grad_freq[1]],[0,0,1]])
   D=np.diag([1./std[0],1./std[1],1./sigma_freq])
   RD=R.dot(D)
   P=RD.dot(RD.T)
   mu=np.array([pos[0],pos[1],freq])
   return (mu,P)

# Create a gassian flux given a WCS with compact support 
def create_gauss_flux(cube,mu,P,peak,cutoff):
   Sigma=P.inv()
   window=np.sqrt(2*np.log(peak/cutoff)*np.diag())
   lower,upper=cube.index_from_window(mu,window)
   feat=cube.get_features(lower,upper)
   C=np.empty_like(feat)
   C[0]=feat[0] - mu[0]
   C[1]=feat[1] - mu[1]
   C[2]=feat[2] - mu[2]
   V=C*(P.dot(C))
   quad=V.sum(axis=0)
   res=np.exp(-quad/2.0)
   res=res.reshape(upper[0]-lower[0],upper[1]-lower[1],upper[2]-lower[2])
   return res,lower,upper

