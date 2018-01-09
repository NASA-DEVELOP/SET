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

pflag = "verbose"
filein = "C:/VIIRS_processing/Clipped Rasters.gdb/VIIRS_2014_09"
myRaster = arcpy.Raster(filein)
imagearr = arcpy.RasterToNumPyArray(myRaster, nodata_to_value = numpy.NaN)
kerneltiffpathhigh = "C:/outputkerneltiffs/kernelubreak30_highlat2.tif"
propagation_arrayhigh = arcpy.RasterToNumPyArray(arcpy.Raster(kerneltiffpathhigh), nodata_to_value = numpy.NaN)
kerneltiffpathlow = "C:/outputkerneltiffs/kernelubreak30_lowlat2.tif"
propagation_arraylow = arcpy.RasterToNumPyArray(arcpy.Raster(kerneltiffpathlow), nodata_to_value = numpy.NaN)
varrprint(propagation_arrayhigh, "High Lat", pflag)
print propagation_arrayhigh[863//2,1259//2]
varrprint(propagation_arraylow, "Low Lat", pflag)
print propagation_arraylow[863//2,1259//2]
propagation_arrayhigh = float32(nan_to_num(propagation_arrayhigh))
propagation_arraylow = float32(nan_to_num(propagation_arraylow))

pad_left = (imagearr.shape[0] - propagation_arrayhigh.shape[0])//2
pad_right = (imagearr.shape[0] - propagation_arrayhigh.shape[0])//2 + 1
pad_up = (imagearr.shape[1] - propagation_arrayhigh.shape[1])//2
pad_down = (imagearr.shape[1] - propagation_arrayhigh.shape[1])//2

pad_leftl = (imagearr.shape[0] - propagation_arraylow.shape[0])//2
pad_rightl = (imagearr.shape[0] - propagation_arraylow.shape[0])//2 + 1
pad_upl = (imagearr.shape[1] - propagation_arraylow.shape[1])//2
pad_downl= (imagearr.shape[1] - propagation_arraylow.shape[1])//2

padded_prophigh = pad(propagation_arrayhigh,((pad_left,pad_right),(pad_up,pad_down)), 'constant', constant_values = 0)
padded_proplow = pad(propagation_arraylow, ((pad_leftl,pad_rightl),(pad_upl,pad_downl)), 'constant', constant_values = 0)
varrprint(padded_prophigh, 'pprop high', pflag)
varrprint(padded_proplow, 'pprop low', pflag)
np_dft_prop_imh = fft.fft2(padded_prophigh)
np_dft_kernel_shifth = fft.fftshift(np_dft_prop_imh)
np_magnitude_spectrumh = 20*log(abs(np_dft_kernel_shifth))

np_dft_prop_iml = fft.fft2(padded_proplow)
np_dft_kernel_shiftl = fft.fftshift(np_dft_prop_iml)
np_magnitude_spectruml = 20*log(abs(np_dft_kernel_shiftl))


plt.subplot(121),plt.imshow(np_magnitude_spectrumh, cmap = 'gray')
plt.title('high conv (size/patterns match)'), plt.xticks([]), plt.yticks([])
plt.subplot(122),plt.imshow(np_magnitude_spectruml, cmap = 'gray')
plt.title('low conv'), plt.xticks([]), plt.yticks([])
plt.show()

filein = "C:/VIIRS_processing/Clipped Rasters.gdb/VIIRS_2014_09"
myRaster = arcpy.Raster(filein)
imagearr = arcpy.RasterToNumPyArray(myRaster, nodata_to_value = numpy.NaN)
viirs_scaling_factor = 10**9
imagearr *= viirs_scaling_factor

np_dft_viirs_im = fft.fft2(imagearr)
np_dft_viirs_shift = fft.fftshift(np_dft_viirs_im)
np_magnitude_spectrum_viirs = 20*log(abs(np_dft_viirs_shift))

kernel_inv_shifth = fft.ifftshift(np_dft_kernel_shifth)
kernel_inv_shiftl = fft.ifftshift(np_dft_kernel_shiftl)
viirs_inv_shift = fft.ifftshift(np_dft_viirs_shift)
varrprint(viirs_inv_shift, 'viirs inv shift', pflag)

FFT_product_inverseh = abs(fft.fftshift(fft.ifft2(kernel_inv_shifth * viirs_inv_shift)))
FFT_product_inversel = abs(fft.fftshift(fft.ifft2(kernel_inv_shiftl * viirs_inv_shift)))
# varrprint(FFT_product_inverseh, 'FFT_product_inverse high', pflag)
# varrprint(FFT_product_inversel, 'FFT_product_inverse low', pflag)
plt.subplot(121),plt.imshow(FFT_product_inverseh, norm = colors.LogNorm(), cmap = 'gray')
plt.title('Conv at High Lat'), plt.xticks([]), plt.yticks([])
plt.subplot(122),plt.imshow(FFT_product_inversel, norm = colors.LogNorm(), cmap = 'gray')
plt.title('Conv at Low Lat'), plt.xticks([]), plt.yticks([])
plt.show()

diff=FFT_product_inverseh - FFT_product_inversel
diffrelative = diff/FFT_product_inverseh
varrprint(FFT_product_inverseh, 'High Lat Convolution', pflag)
varrprint(FFT_product_inversel, 'Low Lat Convolution', pflag)
varrprint(diff, "Difference High Lat - Low Lat", pflag)
varrprint(diffrelative, "Difference normalized to High Lat", pflag)
plt.subplot(121),plt.imshow(diff, norm = colors.LogNorm(), cmap = 'jet')
plt.title('high lat conv - low lat conv'), plt.xticks([]), plt.yticks([])
plt.colorbar()
plt.subplot(122),plt.imshow(diffrelative, norm = colors.LogNorm(), cmap = 'jet')
plt.title('Difference/high lat conv'), plt.xticks([]), plt.yticks([])
plt.colorbar()
plt.show()

