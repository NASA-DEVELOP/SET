import numpy
from numpy import *
from osgeo import gdal, osr

# Read the input raster into a Numpy array
infile = "C:/outputkerneltiffs/kernelubreak30_highlat2.tif"
kerneldata   = gdal.Open(infile)
arr    = kerneldata.ReadAsArray()
# Read VIIRS Data for transform information
imdata = gdal.Open("C:/MonthlyViirs2015/20140901_20140930_75N180W_C.tif")

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


