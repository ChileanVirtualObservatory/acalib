import unittest
import sys
import numpy as np

from astropy.io import fits
from subprocess import call
import shutil
import os

sys.path.append("..")
import acalib.algorithms as acaalgo

from acalib.io import loadFITS_PrimaryOnly
import acalib


"""
	Function to download a fits and load as np.ndarray
"""
def download_and_load():
	if not os.path.exists("files"):
		print("Creating files directory\n#######################")
		os.makedirs("files")

	files = os.listdir("files")
	if not "Orion.methanol.cbc.contsub.image.fits" in files:
		print("Download Fits\n#######################")
		call(["wget", "https://github.com/ChileanVirtualObservatory/bindata/raw/master/fits/cubes/Orion.methanol.cbc.contsub.image.fits"])
		shutil.move("Orion.methanol.cbc.contsub.image.fits","files/")

	hdulist = fits.open("files/Orion.methanol.cbc.contsub.image.fits")
	data = hdulist[0].data
	hdulist.close()
	return data[0,:,:,:]

class TestGMS(unittest.TestCase):
	gms = acaalgo.GMS()
	
	data = download_and_load()
	data = np.sum(data,axis=0)

	def test_nddata(self):
		assert(len(self.gms.run(self.data)) == 1)

class TestCF(unittest.TestCase):
	cf = acaalgo.ClumpFind()

	data = download_and_load()
	data2d = np.sum(data,axis=0)

	def test_nddata_3d(self):
		caa = self.cf.run(self.data)[0]
		assert(caa.min() == 0)
		assert(caa.max() == 18)

	def test_nddata_2d(self):
		caa = self.cf.run(self.data2d)[0]
		assert(caa.min() == 0)
		assert(caa.max() == 3)		


class TestFW(unittest.TestCase):
	fw = acaalgo.FellWalker()

	data = download_and_load()
	data2d = np.sum(data,axis=0)

	def test_nddata_3d(self):
		caa = self.fw.run(self.data)[0]
		assert(caa.min() == 0)
		assert(caa.max() == 10 )

	def test_nddata_2d(self):
		caa = self.fw.run(self.data2d)[0]
		assert(caa.min() == 0)
		assert(caa.max() == 4)


class TestIndexing(unittest.TestCase):
	idx = acaalgo.Indexing({"RANDOM_STATE":1234})

	data = download_and_load()

	#TODO make a better unit test
	def test_run(self):
		result = self.idx.run(self.data).images
		np.testing.assert_equal(len(result),6)


class TestStacking(unittest.TestCase):
	st = acaalgo.Stacking()

	template = download_and_load()
	template = np.sum(template,axis=0)

	cont = acalib.Container()
	cont.images.append(template)
	cont.images.append(template)
	cont.primary = cont.images[0]


	#TODO make a better unit test
	def test_run(self):
		result = self.st.run(self.template,self.cont)
		np.testing.assert_equal(result.shape,(99,99))

if __name__ == '__main__':
    unittest.main()
