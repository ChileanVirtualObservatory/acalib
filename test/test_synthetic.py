import matplotlib.pyplot as plt
import synthetic.imc as imc
import astropy.units as u
import synthetic.vu as vu
import numpy as np
import math

univ=vu.Universe()
univ.create_source('example',0.0*u.deg,1.0*u.deg)
univ.create_source('example2',0.05*u.deg,0.95*u.deg)
model=imc.IMC('Gaussian',np.array([100,70])*u.arcsec,math.pi/9.*u.rad,'33SO2',300*u.deg,10*u.km/u.s,np.array([5,-1])*u.km/(u.s*u.arcsec))
model2=imc.IMC('Gaussian',np.array([30,100])*u.arcsec,math.pi/3.*u.rad,'33SO2',280*u.deg,20*u.km/u.s,np.array([-2,0])*u.km/(u.s*u.arcsec))
model.set_velocity(150*u.km/u.s)
model2.set_velocity(100*u.km/u.s)
univ.add_component('example',model)
univ.add_component('example2',model2)
cube=univ.gen_cube(np.array([0.0,1.0])*u.deg,np.array([1.0,1.0])*u.arcsec, np.array([200,200])*u.arcsec,300*u.GHz, 0.1*u.GHz, 2*u.GHz, 0.0001)
#univ.save_cube(cube,'p33SO2-obs1.fits')
plt.plot(cube.get_spectrum(0.0,1.0))
plt.show()
cube.animate(10,True)

