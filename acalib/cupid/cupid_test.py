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


#config = {'FWHMBEAM':2.0, 'VELORES':2.0, 'ALLOWEDGE':0, 'NAXIS':2, 'RMS':0, 'IDLALG':0, 'MINPIX':10}
config = dict()

"""
CUPID's clumpfind tests
"""
print("--------------------------------------------")
print("PERFORMING CLUMPFIND TESTS")
cf2d_res = clumpfind(data2D, config, rms2D)
cf3d_res = clumpfind(data3D, config, rms3D)
print("Results on 2D data:")
print(type(cf2d_res.data))
print("number of clumps detected: {0}\n".format(np.max(cf2d_res)+1))
for clump in range(np.max(cf2d_res)+1):
    print("number of pixels of clump {0}: {1}".format(clump, np.sum(cf2d_res==clump) ))
print("\nResults on 3D data")
print("number of clumps detected: {0}\n".format(np.max(cf3d_res)+1+1))
for clump in range(np.max(cf3d_res)+1):
    print("number of pixels of clump {0}: {1}".format(clump, np.sum(cf3d_res==clump) ))
print("--------------------------------------------")


"""
CUPID's fellwalker tests
"""
#fw2d = fellwalker(data2D, dict(), rms2D)
#ret2 = fellwalker(data3D, dict(), rms3D)
# ret = []
# for i in range(100):
#     ret.append(clumpfind(data2D, dict(), rms2D).max())
# print ret
