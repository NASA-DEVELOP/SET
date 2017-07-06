"""
Source Name:	Skyglow Estimation Toolbox.pyt
Author:			Wyoming Cross-Cutting II, NASA DEVELOP National Program
Description:	Python toolbox to generate artificial skyglow maps using data from NASA
				and NOAA's Visible Infrared Imaging Radiometer Suite (VIIRS) sensor.

REFERENCES
(1) Falchi, F., P. Cinzano, D. Duriscoe, C.C.M. Kyba, C.D. Elvidge, K. Baugh, B.A. 
	  Portnov, N.A. Rybnikova and R. Furgoni, 2016. The new workd atlas of artificial 
	  night sky brightness. Sci. Adv. 2.
(2) Cinzano, P., F. Falchi, C.D. Elvidge and  K.E. Baugh, 2000. The artificial night sky 
	  brightness mapped from DMSP satellite Operational Linescan System measurements. 
(3) Garstang, R.H., 1989. Night-sky brightness at observatories and sites. Pub. 
	  Astron. Soc. Pac. 10 Mon. Not. R. Astron. Soc. 318.1.
"""

from arcpy import *
from numpy import *
import os.path
import time
import itertools

# Function to streamline parameter creation.
def parameter(dName,name,datatype,paramtype=None,direction=None):
	if paramtype == None:
		paramtype = "Optional"
	if direction == None:
		direction = "Input"
	param = arcpy.Parameter(
		displayName = dName,
		name = name,
		datatype = datatype,
		parameterType = paramtype,
		direction = direction)
	return param

class Toolbox(object):
	# Define the toolbox.
    def __init__(self):
        self.label = "Skyglow Estimation Toolbox"
        self.alias = "skyglow"

        # List of tool classes associated with toolbox.
        self.tools = [CreateArtificialSkyglowMap]

class CreateArtificialSkyglowMap(object):
	# Define the tool.
	def __init__(self):
		self.label = "Create Artificial Skyglow Map"
		self.description = """The Create Artificial Skyglow Map tool receives a kernel file (or produces
							  one based on user inputs) and creates an image of nighttime sky brightness
							  through convolution."""
		self.canRunInBackground = True

	# Define parameters.
	def getParameterInfo(self):

		# Input VIIRS reference
		vImage = parameter("Input VIIRS Reference","vImage","DEFile")

		# Region Latitude (deg), Grand Teton National Park = 43.7904 degrees N
		lat_in = parameter("Region Latitude (deg)","lat","GPDouble")

		# Distance from site to point in sky, relaxing increment (km), u, REF 2, Fig. 6, p. 648
		ubr_in = parameter("Distance from site to night sky (km), relaxing increment","ubr","GPDouble")

		# Zenith angle (deg), z, REF 2, Fig. 6, p. 648
		zen_in = parameter("Zenith Angle (deg)","zen","GPDouble")

		# Azimuth angle (deg)
		azi_in = parameter("Azimuth Angle (deg)","azi","GPDouble")

		# Kernel Output Folder
		kerneloutput = parameter("Output Kernel Folder","kerneloutput","DEFolder",direction="Output")

		# Kernel TIFF File Path
		kerneltiffpath = parameter("Input Kernel File (.tif)","kerneltiffpath","DERasterDataset")

		params = [vImage,lat_in,ubr_in,zen_in,azi_in,kerneloutput,kerneltiffpath]
		return params

	# Source code
	def execute(self,params,messages):
		lat_arg = float(params[1].valueAsText)
		ubr_arg = float(params[2].valueAsText)
		zen_arg = float(params[3].valueAsText)
		azi_arg = float(params[4].valueAsText)
		kerneloutput = params[5].valueAsText
		kerneltiffpath = params[6].valueAsText

		if kerneltiffpath is None:
			# Estimate the 2D propagation function
			propkernel, time = fsum_2d(lat_arg,ubr_arg,zen_arg,azi_arg)
		return

# Function for 2D light propagation
def fsum_2d(lat,ubr,zen,azi):
	AddMessage("**INPUTS**")

	# Central latitude and kernel pixel size
	cent_lat_deg = lat
	cent_lat = cent_lat_deg*pi/180
	pix_deg = .0041666667
	pix = pix_deg*pi/180

	# z, Zenith Angle, REF 2, Fig. 6, p. 648
	AddMessage("Site zenith, z: {} deg".format(zen))

	# Length of u, relaxing integration increment, REF 2, Fig. 6, p. 648
	AddMessage("Length of u for relaxing integration increment, ubr: {} km".format(ubr))

	# Center latitude of region
	AddMessage("Center latitude of region, lat: {} deg".format(cent_lat_deg))

	# Gaussian Earth radius of curvature, REF Wikipedia (https://en.wikipedia.org/wiki/Earth_radius)
	R_T = gauss_earth_curvature_radius(cent_lat)
	AddMessage("Radius of curvature of the Earth calculated from latitude, R_T: {} km".format(R_T))

	# Azimuth angle
	AddMessage("Azimuth angle, azi: {} deg".format(azi))
	AddMessage("******************************")

	# Create latitude/longitude arrays
	source_lat, rltv_long, cent_row, cent_col = create_latlon_arrays(R_T, cent_lat, pix)

	# Distance from source (C) to observation site (O) along ellipsoid surface, REF 2, Fig. 6, p. 648
	# using haversine formula
	D_OC = 2*R_T*arcsin(sqrt(sin((source_lat - cent_lat)/2)**2 + cos(cent_lat)*cos(source_lat)*sin(rltv_long/2)**2))

	# Check of D_OC shape after assigning NaNs outside of 200 km
	AddMessage("Kernel height in pixels, trimmed: {}".format(D_OC.shape[0]))
	AddMessage("Kernel width in pixels, trimmed: {}".format(D_OC.shape[1]))
	widthcenter = (D_OC.shape[1] + 1)//2
	heightcenter = (D_OC.shape[0] +1)//2
	################################## reassignment of center value, need to use better method
	D_OC[cent_row,cent_col] = .01

	# Earth angle from source to site, REF 3, p. 308
	chi = D_OC/R_T

	# beta array, Azimuth angle from line of sight to scatter from site, REF 2, Fig. 6, p. 648
	# http://www.codeguru.com/cpp/cpp/algorithms/article.php/c5115/Geographic-Distance-and-Azimuth-Calculations.htm
	AddMessage("******************************")
	beta = arcsin(sin(pi/2-source_lat)*sin(rltv_long)/sin(chi))
	varrprint(beta, "Beta, Relative azimuth line-of-sight to scatter (rad)")
	test = beta*180/pi
	AddMessage(test[cent_row-2:cent_row+3, cent_col-2:cent_col+3])
	abeta = beta - azi

	# u0, shortest scattering distance based on curvature of the Earth (km), REF 2, Eq. 21, p. 647
	u0 = 2*R_T*sin(chi/2)**2/(sin(zen)*cos(abeta)*sin(chi)+cos(zen)*cos(chi))
	varrprint(u0,"u0")

	# l, Direct line of sight distance between source and observations site (km), REF 2, Appendix A (A1), p. 656
	# l_OC and D_OC are similar as expected
	l_OC = sqrt(4*R_T**2*sin(chi/2)**2)

	# q1, Intermediate quantity (km), REF 2, Appendix A (A1), p. 656, **WITH CORRECTION FROM REF 3, eq. 6, p. 308**
	q1 = R_T*(sin(chi)*sin(zen)*cos(abeta) + cos(chi)*cos(zen) - cos(zen))

	# theta, elevation angle of scatter above source from site (QOC) in radians, REF 2, Appendix A (A1), p. 656
	theta = arccos(q1/l_OC)

	# Get left arrays to cut processing time in half
	abetaleft = abeta[0:,0:widthcenter]
	chileft = chi[0:,0:widthcenter]
	u0left = u0[0:,0:widthcenter]
	l_OCleft = l_OC[0:,0:widthcenter]
	thetaleft = theta[0:,0:widthcenter]

	# Container for propagation array
	PropSumArrayleft = zeros_like(l_OCleft)

	AddMessage("Time for iterations, no threads")
	start = time.time()

	# 2d iteration for integrating from u0 to infinity to create propagation function for each element
	for p,c,u,l,t,b in itertools.izip(nditer(PropSumArrayleft, op_flags=['readwrite']),
									  nditer(chileft, op_flags=['readwrite']),
									  nditer(u0left, op_flags=['readwrite']),
									  nditer(l_OCleft, op_flags=['readwrite']),
									  nditer(thetaleft, op_flags=['readwrite']),
									  nditer(abetaleft, op_flags=['readwrite'])):
		p[...] = fsum_single(R_T, c, u, l, t, b, zen, ubr)
	end = time.time()
	time_sec = end-start
	PropSumArrayright = fliplr(PropSumArrayleft[:,1:])

	# Complete 2d propagation function
	PropSumArray = hstack((PropSumArrayleft, PropSumArrayright))
	return PropSumArray, time_sec

# Calculates the Gaussian Earth radisu of curvature from latitude.
def gauss_earth_curvature_radius(cent_lat_arg):
	# Earth ellipse semi-major orbit, REF Wikipedia (https://en.wikipedia.org/wiki/Earth_radius)
	R_eq = 6378.1370 # km 

	# Earth ellipse semi-minor axis, REF Wikipedia (https://en.wikipedia.org/wiki/Earth_radius)
	R_pol = 6356.7523142 # km

	# Gaussian Earth radius of curvature, REF Wikipedia (https://en.wikipedia.org/wiki/Earth_radius)
	R_curve = ((R_eq**2)*R_pol)/((R_eq*cos(cent_lat_arg))**2 + (R_pol*sin(cent_lat_arg))**2)

	return R_curve

# Function for initial array sizing and create latitude/longitude arrays
def create_latlon_arrays(R_curve, cent_lat, pix):
	p_h = R_curve*pix
	p_w = cos(cent_lat)*R_curve*pix

	# Calculate dimensions for kernel
	kernel_row = int(round(400/p_h))
	kernel_col = int(round(400/p_w))
	if kernel_row%2 == 0:
		kernel_row += 1
	if kernel_col%2 == 0:
		kernel_col += 1
	AddMessage("Kernel height in rows, untrimmed: {}".format(kernel_row))
	AddMessage("Kernel width in columns, untrimmed: {}".format(kernel_col))

	# Create vectors of column and row counts
	col_count = array(range(kernel_col))
	row_count = array(range(kernel_row))

	# Find center column, center row
	center_col = int((kernel_col - 1)/2)
	center_row = int((kernel_row - 1)/2)

	# Create vectors of relative longitude and latitudes (in radians)
	rel_long_vec = array(pix*(col_count - center_col))

	rel_lat_vec = -pix*(row_count - center_row)

	rel_long = tile(rel_long_vec,(kernel_row,1))
	rel_lat = transpose(tile(rel_lat_vec,(kernel_col,1)))
	varrprint(rel_long,'rel_long')
	varrprint(rel_lat,'rel_lat')

	src_lat = rel_lat + cent_lat
	varrprint(src_lat,'src_lat')

	return src_lat, rel_long, center_row, center_col

# For every single pixel in array, takes elements (D_OC, chi, etc.) as arguments
def fsum_single():

def varrprint(varrval, varrtext):
	AddMessage("******************** {} :".format(varrtext))
	AddMessage("{} {}".format(varrval.shape,varrval.dtype))
	AddMessage("{} maximum: {}".format(varrtext, ma.maximum(varrval[~isnan(varrval)])))
	AddMessage("{} minimum: {}".format(varrtext, ma.minimum(varrval[~isnan(varrval)])))