import unittest
from math import pi
from numpy import array, tile, transpose, array_equal, arcsin, sin, cos, sqrt

from darkskypy import (gauss_earth_curvature_radius, sgmapper, create_latlon_arrays,
		       create_beta)


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


class TestBeta(unittest.TestCase):
    def test_some_input(self):
        lat, lon, center_lat, center_lon = create_latlon_arrays(100,0,1)
        d = 2.0*400*arcsin(sqrt(sin((lat - 0)/2.0)**2.0 + cos(0)*cos(lat)*sin(lon/2.0)**2.0))
        d[2,2] = 0.01
        chi = d/400
        correct = array([[2.74724455, 2.77399189, 3.14159265, 3.50919342, 3.53594076],
       			 [3.67005487, 3.63695994, 3.14159265, 2.64622536, 2.61313044],
       			 [4.71238898, 4.71238898, 0, 1.57079633, 1.57079633],
       			 [5.75472309, 5.78781802, 0, 0.49536729, 0.52846221],
       			 [6.67753341, 6.65078607, -0., -0.36760077, -0.39434811]])
	result = create_beta(lat,lon,chi,0,center_lat, center_lon).round(8)
        self.assertTrue(array_equal(result, correct))


if __name__ == '__main__':
    unittest.main()
