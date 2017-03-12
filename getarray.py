# need / to represent true division
from __future__ import division
#single case array build
#import arcpy
import numpy
from numpy import *

#arbitrary, colter bay
central_lat = 43.9047
R_teton = 6367.941
p_deg = .0041666667

#these are in km. They don't seem to quite match up with the 750 meter res of VIIRS...
p_h = R_teton*p_deg*pi/180
p_w = cos(central_lat*pi/180)*R_teton*p_deg*pi/180

#pixel units

kernel_w = round(400/p_w)

kernel_h = round(400/p_h)

if kernel_w%2 == 0:
	kernel_w += 1 
if kernel_h%2 == 0:
	kernel_h += 1

rel_long = zeros((int(kernel_w), int(kernel_h)))
rel_lat = ones_like(rel_long)*central_lat
centralcolumnindex_long = int(round(rel_long.shape[1]/2))-1
centralrowindex_lat = int(round(rel_lat.shape[0]/2))-1

#slice arrays down the middle, then stitch back together after adding element wise
rel_long_left = rel_long.T[0:centralcolumnindex_long]
rel_long_right = rel_long.T[centralcolumnindex_long+1:]

p_h_incr = p_h
for column in rel_long_left:
	column -= p_h_incr
	p_h_incr += p_h
rel_long_left = rel_long_left[ : : -1]
p_h_incr = p_h
for column in rel_long_right:
	column += p_h_incr
	p_h_incr += p_h

centercolumn = rel_long.T[centralcolumnindex_long:centralcolumnindex_long]

rel_long = hstack((rel_long_left.T, centercolumn.T, rel_long_right.T))

# p_w_incr = p_w
# for row in rel_lat:
# 	row += p_w_incr
# 	p_w_incr += p_w

#matrices look right. Columns and rows go to about 400km. Offsets are due to kernel_w and kernel_h rounding
# print rel_long
# print rel_lat

#iterate through and set to distance all distances under 200km radius distance

# distance = 2*R_teton*arcsin(sqrt(sin((0-rel_lat)/2)**2 + cos(rel_lat)*cos(0)*sin((rel_long-0)/2)**2))

# print distance

# nans_dist = zeros_like(rel_long)
# nans_dist[:] = numpy.NAN
# i=0




# in_raster = arcpy.Raster("C:/VIIRS_processing/Clipped Rasters.gdb/VIIRS_2014_06.tif")
# lowerLeft = arcpy.Point(in_raster.extent.XMin, in_raster.extent.YMin)


# haversine_D = 2*R_teton*arcsin(sqrt(sin((lat_O-lat_C)/2)**2 + cos(lat_C)*cos(lat_O)*sin((0)/2)**2))

# arr = arcpy.RasterToNumPyArray(in_raster, nodata_to_value = 0)
# print(arr)
# dtype4array = float
# columns = 1915
# rows = 1428
# K_size = array()
# cellsizecheck = (.0041666667,.0041666667)

