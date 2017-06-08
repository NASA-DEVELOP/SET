import numpy
from numpy import *
from osgeo import gdal, osr
from matplotlib import pyplot as plt
import matplotlib.colors as colors
# Read the input raster into a Numpy array
infile = "./kernelubreak30_highlat2.tif"
kerneldata   = gdal.Open(infile)
arr    = kerneldata.ReadAsArray()
# Read VIIRS Data for transform information
imdata = gdal.Open("./20140901_20140930_75N180W_C.tif")
imarray = imdata.ReadAsArray()

def compare_arr(arr1, arr2, title1, title2):
	plt.subplot(121),plt.imshow(arr1, norm = colors.LogNorm(), cmap = 'gray')
	plt.title(title1), plt.xticks([]), plt.yticks([])
	plt.subplot(122),plt.imshow(arr2, norm = colors.LogNorm(), cmap = 'gray')
	plt.title(title2), plt.xticks([]), plt.yticks([])
	plt.show()

compare_arr(arr, imarray, "Light Scattering Coefficients", "Nighttime Lights VIIRS Image")
# Save out to a GeoTiff

# First of all, gather some information from VIIRS image
[cols,rows] = arr.shape
trans       = imdata.GetGeoTransform()
proj        = imdata.GetProjection()
nodatav     = 0
outfile     = 'C:/outputkerneltiffs/gdalproduct_referenced.tif'

# Create the georeffed file, using the information from the VIIRS image
outdriver = gdal.GetDriverByName("GTiff")
outdata   = outdriver.Create(str(outfile), rows, cols, 1, gdal.GDT_Float32)

# Write data to the file, which is the kernel array in this example
outdata.GetRasterBand(1).WriteArray(arr)

# Set a no data value
outdata.GetRasterBand(1).SetNoDataValue(nodatav)

# Georeference the image
outdata.SetGeoTransform(trans)

# Write projection information
outdata.SetProjection(proj)


