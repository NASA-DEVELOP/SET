# need / to represent true division
from __future__ import division
#single case array build
#import arcpy
import numpy
from numpy import *

#arbitrary radius and center latitude, for testing purposes
R_teton = 6367.941
cent_lat = 40.0
p_deg = .0041666667

def get_dist_arr(central_lat, radius, pixel_deg):
	#these are in km. They don't seem to quite match up with the 750 meter res of VIIRS...
	p_h = radius*pixel_deg*pi/180
	p_w = cos(central_lat*pi/180)*radius*pixel_deg*pi/180
	central_lat_km = central_lat*radius*pi/180
	
	#pixel units
	kernel_w = round(400/p_w)
	print "kernel width in pixels: {}".format(kernel_w)
	kernel_h = round(400/p_h)
	print "kernel height in pixels: {}".format(kernel_h)
	if kernel_w%2 == 0:
		kernel_w += 1 
	if kernel_h%2 == 0:
		kernel_h += 1

	#made relative arrays the same size, easier, few extra distance calculations will be trimmed
	rel_long = zeros((int(kernel_w), int(kernel_w)))
	rel_lat = ones((int(kernel_w), int(kernel_w)))*central_lat_km
	

	centralcolumnindex_long = int(round(rel_long.shape[1]/2))-1
	centralrowindex_lat = int(round(rel_lat.shape[0]/2))-1

	#slice arrays down the middle, then stitch back together after adding pixel increment element wise
	def slice_n_stack(arr, pix_measure, central_index):
		arr_left = arr.T[0:central_index]
		arr_right = arr.T[central_index+1:]
		
		p_incr = pix_measure
		for column in arr_left:
			column -= p_incr
			p_incr += pix_measure
		arr_left = arr_left[ : : -1]
		
		p_incr = pix_measure
		for column in arr_right:
			column += p_incr
			p_incr += pix_measure
		

		centercolumn = arr.T[central_index:central_index+1]
		
		result = hstack((arr_left.T, centercolumn.T, arr_right.T))
		return result
		

	rel_long = slice_n_stack(rel_long, p_h, centralcolumnindex_long)
	rel_lat = slice_n_stack(rel_lat.T, p_w, centralrowindex_lat)
	
	# convert to radians for haversine formula
	rel_long = rel_long/radius
	rel_lat = rel_lat/radius
	central_lat = central_lat*pi/180
	distance = 2*radius*arcsin(sqrt(sin((central_lat-rel_lat)/2)**2 + cos(rel_lat)*cos(central_lat)*sin((rel_long-0)/2)**2))
	#wooohoo makes sense! longest distance should be about 282 from center to hypotenuse.
	return distance


distance_kernel = get_dist_arr(cent_lat, R_teton, p_deg)
for x in nditer(distance_kernel, op_flags=['readwrite']):
	if x > 201:
		x[...] = numpy.NaN
		
distance_kernel = distance_kernel[:, ~isnan(distance_kernel).any(axis=0)]
Chi_kernel = distance_kernel/R_teton


# in_raster = arcpy.Raster("C:/VIIRS_processing/Clipped Rasters.gdb/VIIRS_2014_06.tif")
# lowerLeft = arcpy.Point(in_raster.extent.XMin, in_raster.extent.YMin)

# arr = arcpy.RasterToNumPyArray(in_raster, nodata_to_value = 0)
# print(arr)
# dtype4array = float
# columns = 1915
# rows = 1428
# K_size = array()
# cellsizecheck = (.0041666667,.0041666667)

