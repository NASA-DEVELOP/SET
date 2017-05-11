from __future__ import division
import numpy
from numpy import *
import itertools
import time
import archook
archook.get_arcpy()
import arcpy
import threading
from scipy import ndimage
import os.path
from matplotlib import pyplot as plt
import matplotlib.colors as colors
def varrprint(varrval, varrtext, print_flag):
	if print_flag != 'quiet':
		print '********************',
		print varrtext,
		print ':'
		print varrval.shape,
		print varrval.dtype
		print varrtext, 
		print 'maximum: {}'.format(ma.maximum(varrval[~isnan(varrval)])),
		print ',',
		print varrtext,
		print 'minimum: {}'.format(ma.minimum(varrval[~isnan(varrval)]))

pflag = "verbose"
kerneltiffpath = "C:/outputkerneltiffs/kernelubreak30.tif"
propagation_array1 = arcpy.RasterToNumPyArray(arcpy.Raster(kerneltiffpath), nodata_to_value = numpy.NaN)
filein = "C:/VIIRS_processing/Clipped Rasters.gdb/VIIRS_2014_08"
myRaster = arcpy.Raster(filein)
imagearr = arcpy.RasterToNumPyArray(myRaster, nodata_to_value = numpy.NaN)
varrprint(imagearr, 'VIIRS', pflag)

viirs_scaling_factor = 10**9
imagearr *= viirs_scaling_factor

propagation_array1 = float32(nan_to_num(propagation_array1))

# need to be same size to multiple fourier transforms
padded_prop = pad(propagation_array1,((282,283),(328,328)), 'constant', constant_values = 0)
varrprint(padded_prop, 'kernel scaled and padded',pflag)
varrprint(imagearr, 'VIIRSscaled', pflag)

np_dft_prop_im = fft.fft2(padded_prop)
np_dft_kernel_shift = fft.fftshift(np_dft_prop_im)
np_magnitude_spectrum = 20*log(abs(np_dft_kernel_shift))
varrprint(np_dft_prop_im, 'np_dft_kernel', pflag)

np_dft_viirs_im = fft.fft2(imagearr)
np_dft_viirs_shift = fft.fftshift(np_dft_viirs_im)
np_magnitude_spectrum_viirs = 20*log(abs(np_dft_viirs_shift))
varrprint(np_dft_viirs_im, 'np_dft_viirs', pflag)

kernel_inv_shift = fft.ifftshift(np_dft_kernel_shift)
viirs_inv_shift = fft.ifftshift(np_dft_viirs_shift)

FFT_product_inverse = abs(fft.fftshift(fft.ifft2(kernel_inv_shift * viirs_inv_shift)))
varrprint(FFT_product_inverse, 'FFT_product_inverse', pflag)

lwr_l_long, lwr_l_lat = -112.617213, 41.417855
lower_left = arcpy.Point(lwr_l_long, lwr_l_lat)
cent_lat = 46.755666
p_deg = .0041666667
x_size = cos(cent_lat*pi/180)*p_deg
y_size = p_deg

output = arcpy.NumPyArrayToRaster(FFT_product_inverse, lower_left, x_size, y_size, 0)
geotiffpath = "C:/outputkerneltiffs/testarcpyproduct.tif"
output.save(geotiffpath)