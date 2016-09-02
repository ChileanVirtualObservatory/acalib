import sys
sys.path.append('../../')
from acalib import *
from acalib.cupid import clumpfind, fellwalker


"""
loading data
"""
binpath = '../../bindata/fits/cubes/'

orion_path = binpath + 'ALMA01000740.fits'
container = load_fits(orion_path)
orion = container.primary


data3D = orion.data.astype(np.float64)
rms3D = np.sqrt((data3D*data3D).sum()/data3D.size)

data2D = data3D.sum(axis=0)
rms2D = np.sqrt((data2D*data2D).sum()/data2D.size)


"""
CUPID's clumpfind call
"""
ret1 = fellwalker(data2D, dict(), rms2D)
ret2 = fellwalker(data3D, dict(), rms3D)
print ret1.max()
print ret2.max()
# ret = []
# for i in range(100):
#     ret.append(clumpfind(data2D, dict(), rms2D).max())
# print ret
