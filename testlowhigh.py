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
kerneltiffpathhigh = "C:/outputkerneltiffs/kernelubreak30_highlat2.tif"
propagation_arrayhigh = arcpy.RasterToNumPyArray(arcpy.Raster(kerneltiffpathhigh), nodata_to_value = numpy.NaN)
kerneltiffpathlow = "C:/outputkerneltiffs/kernelubreak30_lowlat.tif"
propagation_arraylow = arcpy.RasterToNumPyArray(arcpy.Raster(kerneltiffpathlow), nodata_to_value = numpy.NaN)
varrprint(propagation_arrayhigh, "High Lat", pflag)
print propagation_arrayhigh[863//2,1259//2]
varrprint(propagation_arraylow, "Low Lat", pflag)
print propagation_arraylow[863//2,1259//2]

plt.subplot(121),plt.imshow(propagation_arrayhigh, norm = colors.LogNorm(), cmap = 'gray')
plt.title('High Lat'), plt.xticks([]), plt.yticks([])
plt.subplot(122),plt.imshow(propagation_arraylow, norm = colors.LogNorm(), cmap = 'gray')
plt.title('Low Lat'), plt.xticks([]), plt.yticks([])
plt.show()
