import astropy.units as u
import numpy as np
import core.parameter as par

def clump_to_gauss(pos,std,angle,freq,fwhm,gradient,equiv=u.doppler_radio):
   # Parameter sanitization
   pos=par.to_deg(pos)
   std=par.to_deg(std)
   angle=par.to_rad(angle)
   freq=par.to_hz(freq)
   #print "fwhm",fwhm
   #print "freq",freq
   sigma=par.fwhm_to_sigma(freq - par.vel_to_freq(fwhm,freq,equiv))
   #print "sigma",sigma
   grad= freq/u.deg -  par.to_hz_deg(gradient,freq,equiv) 
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

# Create a gassian flux for a cube with compact support 
def create_gauss_flux(cube,mu,P,peak,cutoff):
   Sigma=np.linalg.inv(P)
   window=np.sqrt(2*np.log(peak/cutoff)*np.diag(Sigma))
   #print "win",window
   lower,upper=cube.index_from_window(mu,window)
   print lower,upper
   feat=cube.get_features(lower,upper)
   C=np.empty_like(feat)
   C[0]=feat[0] - mu[0]
   C[1]=feat[1] - mu[1]
   C[2]=feat[2] - mu[2]
   V=C*(C.dot(P))
   quad=V.sum(axis=1)
   res=np.exp(-quad/2.0)*peak
   print res.max()
   res=res.reshape(upper[0]-lower[0],upper[1]-lower[1],upper[2]-lower[2])
   return res,lower,upper

