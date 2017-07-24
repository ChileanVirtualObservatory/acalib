import unittest
import sys
import numpy as np
sys.path.append("../..")
import acalib.core.statistics as acasta


class TestStatistics(unittest.TestCase):
    def test_rms(self):
        np.random.seed(0)
        data2d = np.random.rand(2,2)
        np.testing.assert_almost_equal(acasta.rms(data2d),0.60681823112676958)
        
        np.random.seed(0)
        data3d = np.random.rand(2,2,2)
        np.testing.assert_almost_equal(acasta.rms(data3d),0.6180936123671793)


    def test_snr_estimation(self):
        np.random.seed(0)
        data3d = np.random.rand(2,2,2)
        np.testing.assert_almost_equal(acasta.snr_estimation(data3d),1.046)


if __name__ == '__main__':
    unittest.main()     
