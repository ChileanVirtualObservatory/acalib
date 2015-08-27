import matplotlib.pyplot as plt
import synthetic.imc as imc
import astropy.units as u
import synthetic.vu as vu
import numpy as np
import math

#TODO This is not working for near zero positions!!! (it is fault of the wcs i think!)
univ=vu.Universe()

# Create Source
center=[1.0,1.0]*u.deg
univ.create_source('example',center)

# Defines a new component
mol_list=dict()
mol_list['33SO2']=[1.0,10.0]* u.Jy/u.beam
temp=300*u.K
offset=np.array([0,0])*u.arcsec
std = np.array([10,7])*u.arcsec
angle=math.pi/9.*u.rad
fwhm=30*u.km/u.s
gradient=np.array([0.0,0.0])*u.km/(u.s*u.arcsec)
rad_vel=150*u.km/u.s
# Create Component
model=imc.GaussianIMC(mol_list,temp,offset,std,angle,fwhm,gradient)
model.set_velocity(rad_vel)
univ.add_component('example',model)

# Defines a new component
mol_list=dict()
mol_list['33SO2']=[0.5,5.0]* u.Jy/u.beam
temp=300*u.K
offset=np.array([20,-20])*u.arcsec
std = np.array([4,12])*u.arcsec
angle=math.pi/4.*u.rad
fwhm=10*u.km/u.s
gradient=np.array([-3.0,3.0])*u.km/(u.s*u.arcsec)
rad_vel=150*u.km/u.s
# Create Component
model=imc.GaussianIMC(mol_list,temp,offset,std,angle,fwhm,gradient)
model.set_velocity(rad_vel)
univ.add_component('example',model)

# Defines a new component
mol_list=dict()
mol_list['33SO2']=[0.5,8.0]* u.Jy/u.beam
temp=300*u.K
offset=np.array([-20,20])*u.arcsec
std = np.array([4,12])*u.arcsec
angle=math.pi/6.*u.rad
fwhm=10*u.km/u.s
gradient=np.array([3.0,-3.0])*u.km/(u.s*u.arcsec)
rad_vel=150*u.km/u.s
# Create Component
model=imc.GaussianIMC(mol_list,temp,offset,std,angle,fwhm,gradient)
model.set_velocity(rad_vel)
univ.add_component('example',model)


# Create Cube
ang_res=np.array([1.0,1.0])*u.arcsec
fov=np.array([200,200])*u.arcsec
freq=300*u.GHz
spe_res=0.005*u.GHz
bw=2*u.GHz
noise=0.001*u.Jy/u.beam

(cube,tab)=univ.gen_cube(center,ang_res,fov, freq,spe_res,bw,noise)
print tab
plt.plot(cube.get_stacked(axis=(1,2)))
plt.show()
plt.clf
plt.imshow(cube.get_stacked())
plt.show()

