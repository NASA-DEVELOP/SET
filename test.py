import unittest
from math import pi

from darkskypy import gauss_earth_curvature_radius, sgmapper


class TestGaussianEarthCurvatureRadius(unittest.TestCase):
    def test_0(self):
        self.assertEqual(gauss_earth_curvature_radius(0), 6356.7523142)

    def test_30(self):
        self.assertEqual(gauss_earth_curvature_radius(30*pi/180), 6367.408777700097)


class TestSGMapper(unittest.TestCase):
    def test_OOB_lat(self):
        with self.assertRaises(ValueError):
            sgmapper(5000, 1, 0, 0, 'fake_file')


if __name__ == '__main__':
    unittest.main()
