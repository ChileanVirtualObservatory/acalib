import unittest
import sys
import numpy as np
sys.path.append("../..")
import acalib.core.analysis as acaana


class TestAnalysis(unittest.TestCase):
    data =  np.load("test/data.npy")
    def test_rms(self):
        data2d = self.data[0,:2,:2]
        np.testing.assert_almost_equal(acaana.rms(data2d),0.82766225546)

        data3d = self.data[:2,:2,:2]
        np.testing.assert_almost_equal(acaana.rms(data3d),0.834985136799)


    def test_snr_estimation(self):
        data3d = self.data[:2,:2,:2]
        np.testing.assert_almost_equal(acaana.snr_estimation(data3d),1.072)


    def test_integrate(self):
        data3d = self.data[:2,:2,:2]

        result0axis = np.array([[ 1.42758652,  1.50110189],
                                [ 1.84584797,  1.85519773]])

        np.testing.assert_almost_equal(acaana.integrate(data3d,axis=0),result0axis)

    def test_spectra_sketch(self):
        data_spectra = self.data
        result =(
            np.array([ 0.89777733,  1.40310624,  1.17576868,  2.53032862,  0.,          0. ,         0.,
  0.,          0. ,         0.        ]),
            [slice(0, 4, None)])

        spectra,slices=acaana.spectra_sketch(data_spectra,10,random_state=1)
        np.testing.assert_almost_equal(spectra, result[0])


    def test_morph(self):
        data = self.data[:,0,1]

        acum = np.cumsum(data)
        diff = acaana.differenceImpl(acum)
        result_diff = np.array([0.71920952,  0.78189237,  0.96274427,  1.26236311,  1.07830751,  1.35241923,
  1.62192038,  2.34231765,  1.82180687,  2.41933219])
        np.testing.assert_almost_equal(diff,result_diff)

        boxing = acaana.segmentationImpl(diff)
        result_boxing = np.array([0.,  1.,  1.,  0.,  0.,  1.,  1.,  0.,  0.,  0.])
        np.testing.assert_almost_equal(boxing,result_boxing)

        boxing = acaana.erosionImpl(boxing)
        result_erosion = np.array([ 1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  0.,  0.])
        np.testing.assert_almost_equal(boxing,result_erosion)


        x = np.array([1,2,3,4,5])
        y = np.array([0,1,0,1,0])
        res = np.array([0,2,0,4,0])
        np.testing.assert_almost_equal(acaana._masking(x,y),res)

    def test_pixel_processing(self):
        data = self.data[:,0,1]
        result_processing = np.array([  0.71920952,  0.78189237,  0.24353474,  0.48047074,  0.11556325, 0.09005613,
  0.54361286,  0.98989842,  0.,          0.,        ])
        np.testing.assert_almost_equal(acaana._pixel_processing(data),result_processing)


    def test_optimal_w(self):
        data = self.data[0,:,:]
        np.testing.assert_almost_equal(acaana._optimal_w(data, p=0.3),16)

    def test_vel_stacking(self):
        data = self.data[:,:5,:5]
        result = np.array([[ 1.12021683,  0.59603398,  1.32836238,  0.73973993,  1.53711335],
       [ 1.89774464,  0.24789119,  1.50520932,  0.7577545 ,  0.04810489],
       [ 1.17531154,  1.39070063,  0.24306376,  0.77872804,  1.42740204],
       [ 0.51866241,  0.48976859,  1.7431068 ,  0.21277142,  1.5490394 ],
       [ 1.14947174,  1.57444903,  0.72570653,  1.5528936 ,  0.76760598]])

        np.testing.assert_almost_equal(acaana.vel_stacking(data,slice(3,5)), result)

if __name__ == '__main__':
    unittest.main()
