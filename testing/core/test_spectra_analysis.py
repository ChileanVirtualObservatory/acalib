import unittest
import sys
import numpy as np
sys.path.append("../..")
import acalib.core.spectra_analysis as acaana


class TestSpectraAnalysis(unittest.TestCase):

    def test_spectra_sketch(self):
        random = np.random.RandomState(0)
        data_spectra = random.rand(10,100,100)
        result =(
            np.array([2.194053 ,  4.3325863,  4.220324 ,  2.2198028,  2.9155446,
        0.       ,  0.       ,  0.       ,  0.       ,  0.            ]),
            [slice(0, 9, None)])

        spectra,slices=acaana.spectra_sketch(data_spectra,10,random_state=1)
        np.testing.assert_almost_equal(spectra, result[0])

    def test_peakdet(self):
        random = np.random.RandomState(0)
        data_spectra = random.rand(100)
        result = (np.array([[  8.        ,   0.96366276],
       [ 20.        ,   0.97861834],
       [ 52.        ,   0.98837384]]), np.array([[  1.60000000e+01,   2.02183974e-02],
       [  3.40000000e+01,   1.87898004e-02]]))
        peaks_max,peaks_min=acaana.peakdet(data_spectra,0.9)
        np.testing.assert_almost_equal(peaks_max, result[0])
        np.testing.assert_almost_equal(peaks_min, result[1])


if __name__ == '__main__':
    unittest.main()     
