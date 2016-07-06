import sys
sys.path.append('../../')

import matplotlib.pyplot as plt
import acalib.synthetic.imc as imc
import astropy.units as u
from acalib.synthetic import *
from acalib.io import *
import numpy as np
import math



#TODO This is not working for near zero positions!!! (it is fault of the wcs i think!)
univ=Universe()

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
model=GaussianIMC(mol_list,temp,offset,std,angle,fwhm,gradient)
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
model=GaussianIMC(mol_list,temp,offset,std,angle,fwhm,gradient)
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
model=GaussianIMC(mol_list,temp,offset,std,angle,fwhm,gradient)
model.set_velocity(rad_vel)
univ.add_component('example',model)


# Create Cube
ang_res=np.array([3.0,3.0])*u.arcsec
fov=np.array([200,200])*u.arcsec
freq=299.898*u.GHz
spe_res=0.002*u.GHz
bw=0.2*u.GHz
noise=0.001*u.Jy/u.beam

cont = univ.gen_cube(center, ang_res, fov, freq, spe_res, bw, noise,noise/50.0)


volume3D(cont.primary)
contour3D(cont.primary)

