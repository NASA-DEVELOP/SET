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
from osgeo import gdal, osr
import skimage.external.tifffile

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

# pflag = "verbose"
# kerneltiffpathlow = "C:/outputkerneltiffs/kernelubreak30_lowlat.tif"
# propagation_array1 = arcpy.RasterToNumPyArray(arcpy.Raster(kerneltiffpathlow), nodata_to_value = numpy.NaN)
# varrprint(propagation_array1, 'kernel',pflag)
# filein = "C:/VIIRS_processing/Clipped Rasters.gdb/VIIRS_2014_09"
# myRaster = arcpy.Raster(filein)
# imagearr = arcpy.RasterToNumPyArray(myRaster, nodata_to_value = numpy.NaN)
# varrprint(imagearr, 'VIIRS', pflag)

# viirs_scaling_factor = 10**9
# imagearr *= viirs_scaling_factor

# propagation_array1 = float32(nan_to_num(propagation_array1))

# # need to be same size to multiple fourier transforms
# padded_prop = pad(propagation_array1,((282,283),(328,328)), 'constant', constant_values = 0)
# varrprint(padded_prop, 'kernel scaled and padded',pflag)
# varrprint(imagearr, 'VIIRSscaled', pflag)

# np_dft_prop_im = fft.fft2(padded_prop)
# np_dft_kernel_shift = fft.fftshift(np_dft_prop_im)
# np_magnitude_spectrum = 20*log(abs(np_dft_kernel_shift))
# varrprint(np_dft_prop_im, 'np_dft_kernel', pflag)

# np_dft_viirs_im = fft.fft2(imagearr)
# np_dft_viirs_shift = fft.fftshift(np_dft_viirs_im)
# np_magnitude_spectrum_viirs = 20*log(abs(np_dft_viirs_shift))
# varrprint(np_dft_viirs_im, 'np_dft_viirs', pflag)

# kernel_inv_shift = fft.ifftshift(np_dft_kernel_shift)
# viirs_inv_shift = fft.ifftshift(np_dft_viirs_shift)

# FFT_product_inverse = abs(fft.fftshift(fft.ifft2(kernel_inv_shift * viirs_inv_shift)))
# varrprint(FFT_product_inverse, 'FFT_product_inverse', pflag)

# # lwr_l_long, lwr_l_lat = -114.6687495, 40.8354164
# # lower_left = arcpy.Point(lwr_l_long, lwr_l_lat)
# # cent_lat = 40.8797
# # p_deg = .0041666667
# # x_size = cos(cent_lat*pi/180)*p_deg
# # y_size = p_deg

unreffed_tiff = "C:/outputkerneltiffs/kernelubreak30_lowlat.tif"
# skimage.external.tifffile.imsave(src_filename, FFT_product_inverse)

reffed_tiff = 'C:/outputkerneltiffs/gdalproduct_lowlat_referenced.tif'

# Opens source dataset
src_ds = gdal.Open(unreffed_tiff)
format = "GTiff"
driver = gdal.GetDriverByName(format)

# Open destination dataset
dst_ds = driver.CreateCopy(reffed_tiff, src_ds, 0)

# Specify raster location through geotransform array
# (uperleftx, scalex, skewx, uperlefty, skewy, scaley)
# Scale = size of one pixel in units of raster projection
# this example below assumes 100x100
gt = [-114.6687495, .0041666667, 0, 46.7854164, 0, -.0041666667]

# Set location
dst_ds.SetGeoTransform(gt)

# Get raster projection
epsg = 4326
srs = osr.SpatialReference()
srs.ImportFromEPSG(epsg)
dest_wkt = srs.ExportToWkt()

# Set projection
dst_ds.SetProjection(dest_wkt)

# Close files
dst_ds = None
src_ds = None