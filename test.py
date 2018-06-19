import unittest
from math import pi
from numpy import array, tile, transpose, array_equal

from darkskypy import gauss_earth_curvature_radius, sgmapper, create_latlon_arrays


class TestGaussianEarthCurvatureRadius(unittest.TestCase):
    """Test gauss_earth_curvature_radius functionality on different inputs"""
    def test_0(self):
        self.assertEqual(gauss_earth_curvature_radius(0), 6356.7523142)

    def test_30(self):
        self.assertEqual(gauss_earth_curvature_radius(30*pi/180), 6367.408777700097)


class TestSGMapper(unittest.TestCase):
    def test_OOB_lat(self):
        with self.assertRaises(ValueError):
            sgmapper(5000, 1, 0, 0, 'fake_file')


class TestLatLonArrays(unittest.TestCase):
    def test_some_input(self):
        result = create_latlon_arrays(100,0,1)
        self.assertTrue(array_equal(result[0], transpose(tile((array(range(5))-2)*-1, (5,1)))))
        self.assertTrue(array_equal(result[1], tile(array(range(5))-2, (5,1))))
        self.assertEqual(result[2], 2)
        self.assertEqual(result[3], 2)


if __name__ == '__main__':
    unittest.main()
