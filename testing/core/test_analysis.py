import unittest
import sys
import numpy as np
sys.path.append("../..")
import acalib.core.analysis as acaana


class TestAnalysis(unittest.TestCase):
    def test_rms(self):
        np.random.seed(0)
        data2d = np.random.rand(2,2)
        np.testing.assert_almost_equal(acaana.rms(data2d),0.60681823112676958)
        
        np.random.seed(0)
        data3d = np.random.rand(2,2,2)
        np.testing.assert_almost_equal(acaana.rms(data3d),0.6180936123671793)


    def test_snr_estimation(self):
        np.random.seed(0)
        data3d = np.random.rand(2,2,2)
        np.testing.assert_almost_equal(acaana.snr_estimation(data3d),1.046)


    def test_integrate(self):
        np.random.seed(0)
        data3d = np.random.rand(2,2,2)

        result0axis = np.array([[ 0.9724683 ,  1.36108348],
                                [ 1.04035059,  1.43665618]])

        result1axis = np.array([[ 1.15157688,  1.26007255],
                               [ 0.86124201,  1.53766711]])

        result2axis = np.array([[ 1.26400287,  1.14764656],
                                [ 1.06954891,  1.32936021]])

        np.testing.assert_almost_equal(acaana.integrate(data3d,axis=0),result0axis)
        np.testing.assert_almost_equal(acaana.integrate(data3d,axis=1),result1axis)
        np.testing.assert_almost_equal(acaana.integrate(data3d,axis=2),result2axis)


    def test_spectra_sketch(self):
        random = np.random.RandomState(0)
        data_spectra = random.rand(10,100,100)
        result =(
            np.array([2.194053 ,  4.3325863,  4.220324 ,  2.2198028,  2.9155446,
        0.       ,  0.       ,  0.       ,  0.       ,  0.            ]),
            [slice(0, 9, None)])

        spectra,slices=acaana.spectra_sketch(data_spectra,10,random_state=1)
        np.testing.assert_almost_equal(spectra, result[0])


    def test_morph(self):
        random = np.random.RandomState(0)
        data = random.rand(10)

        
        acum = np.cumsum(data)
        diff = acaana.differenceImpl(acum)

        result_diff = np.array([ 0.5488135 ,  0.71518937,  1.15157688,  1.26007255,  1.57523168,
        1.90596666,  2.01281889,  2.79773966,  2.97648165,  3.18118118])
        np.testing.assert_almost_equal(diff,result_diff)

        boxing = acaana.segmentationImpl(diff)
        result_boxing = np.array([  6.92472057e-310,   1.00000000e+000,   1.00000000e+000,
         1.00000000e+000,   1.00000000e+000,   1.00000000e+000,
         1.00000000e+000,   1.00000000e+000,   1.00000000e+000,
         6.92465131e-310])
        np.testing.assert_almost_equal(boxing,result_boxing)

        boxing = acaana.erosionImpl(boxing)
        result_erosion = np.array([ 1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.])
        np.testing.assert_almost_equal(boxing,result_erosion)


        x = np.array([1,2,3,4,5])
        y = np.array([0,1,0,1,0])
        res = np.array([0,2,0,4,0])
        np.testing.assert_almost_equal(acaana._masking(x,y),res)

    def test_pixel_processing(self):
        random = np.random.RandomState(0)
        
        data = random.rand(10)
        result_processing = np.array([ 0.5488135 ,  0.71518937,  0.60276338,  0.54488318,  0.4236548 ,
        0.64589411,  0.43758721,  0.891773  ,  0.96366276,  0.38344152])
        np.testing.assert_almost_equal(acaana._pixel_processing(data),result_processing)


    def test_optimal_w(self):
        random = np.random.RandomState(0)
        data = random.rand(100,100)
        np.testing.assert_almost_equal(acaana._optimal_w(data, p=0.3),25)

    def test_vel_stacking(self):
        random = np.random.RandomState(0)
        data = random.rand(10,5,5)
        result = np.array([[ 0.71700433,  0.55281494,  0.85539058,  1.25832874,  0.36748086],
       [ 0.89414051,  1.00630493,  0.6363994 ,  0.91555375,  1.51935047],
       [ 0.71251487,  1.36965673,  0.79341979,  0.87338345,  1.74309402],
       [ 0.71507469,  1.54851358,  0.71307074,  1.59806257,  0.98193768],
       [ 0.90844564,  1.08783732,  0.97619118,  1.47293023,  0.42855052]])

        np.testing.assert_almost_equal(acaana.vel_stacking(data,slice(3,5)), result)

if __name__ == '__main__':
    unittest.main()     
