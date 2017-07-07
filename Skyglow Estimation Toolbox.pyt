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
from osgeo import gdal

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
			array_to_geotiff(propkernel, kerneltiffpath)
			varrprint(propkernel,"propagation array")
			AddMessage("Time for propagation function ubreak 10: {}".format(time))
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

	# Test array subsets to reduce processing time.
	chileft = chi[427:432,530:widthcenter]
	u0left = u0[427:432,530:widthcenter]
	l_OCleft = l_OC[427:432,530:widthcenter]
	thetaleft = theta[427:432,530:widthcenter]
	abetaleft = abeta[427:432,530:widthcenter]

	# Container for propagation array
	PropSumArrayleft = zeros_like(l_OCleft)

	AddMessage("Time for iterations, no threads")
	start = time.time()

	##################### VERONICA: PERCENT COMPLETE/PROGRESS BAR CODE

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
def fsum_single(R_T,chi,u0,l_OC,theta,beta_farg,zen_farg,ubrk_farg,K_am_arg = 1.0, du_farg = .2):
	if isnan(l_OC):
		return nan

	# Scattering distance increment, delta u (km)
	du0 = du_farg
	# Zenith angle at observation site (rad)
	zen = zen_farg
	# Horizontal angle at observation site from line-of-sight (OC) to line of incoming scattered light (OQ) (rad)
	beta = beta_farg
	# Length of u, scattering ditance, at which the integration increment is relaxed (km)
	ubrk = ubrk_farg
	# Parameter for relative importance of aerosols to molecules, REF 1, p.10
	K_am = K_am_arg

	# Constants:
	# N_m,0: Molecular density at sea level (cm^-3), REF 2, p.645
	N_m0 = 2.55e19
	# c: Inverse scale altitude (km^-1), REF 2, p.645
	c_isa = 0.104
	# sigma_m: Integrated Rayleigh scattering cross-section in the V-band (cm^-2*sr^-1), REF 2, p. 646
	sig_m = 1.136e-26

	# Falchi best fit parameters for normalized emission
	# Weights for angular distributions, REF 1, p.21
	W_a = 1.9e-3
	W_b = 5.2e-4
	W_c = 7.6e-5
	# Light loss per night hour after midnight
	d_ll = -0.045

	# Containers for variables that update in loop
	u_OQ = u0 # km
	total_sum = 0
	df_prop = 1

	# Total Propagation stable to 3 significant figures
	stability_limit = 0.001

	# Loop counter
	lc = 1

	# START OF "u" LOOP, REF 2, Eq. (3)
	while u_OQ < 30.0:
		if u_OQ < ubrk:
			du = du0
		else:
			du = du0 * 10.0

		# Distance from source to scattering (CQ) (km), s, REF 2, Appendix A (A1), p.656
		# Equation is wrong in REF 2 (Cinzano). Changed angle from chi to theta, REF 3, p. 308, Equation 7 (Garstang)
		s_CQ = sqrt((u_OQ - l_OC)**2 + 4*u_OQ*l_OC*sin(theta/2)**2)

		# Height of scattering above Earth reference surface Q (km), h, REF 2, Appendix A (A1), p. 656
		h_Q = R_T*(sqrt(1 + (u_OQ**2 + 2*u_OQ*R_T*cos(zen))/R_T**2) - 1)

		# Elevation angle of emission from C to Q over line from C to O (QCO) (rad), phi, REF 2, Appendix A (A1), p.656
		phi = arccos((l_OC**2 + s_CQ**2 - u_OQ**2)/(2*l_OC*s_CQ))
		phi_deg = phi*180/pi

		# Scattering angle at Q (rad), omega, REF 2, Appendix A (A1), p. 656
		omega = theta + phi
		omega_deg = omega*180/pi

		# Intermediate quantity (km), q3, REF 2, Appendix A (A1), p. 656
		q3 = u_OQ*cos(zen)*cos(chi) - 2*R_T*sin(chi/2)**2
		# Intermediate quantity (km), q2, (q2=q3 if z=0), REF 2, Appendix A (A1), p. 656
		q2 = u_OQ*sin(zen)*cos(beta)*sin(chi) + q3

		# Emission angle from source (rad), psi, REF 2, Appendix A (A1), p. 656
		psi = arccos(q2/s_CQ)
		psi_deg = psi*180/pi

		# Scale height of aerosols (km^-1), a, REF 2, p. 646
		a_sha = 0.657 + 0.059*K_am

		# Intermediate quantity u-path, p4, REF 2, Appendix A (A2), p. 656
		p4 = (a_sha**2*u_OQ**2*cos(zen)**2 + 2*a_sha*u_OQ*cos(zen) + 2)*exp(-a_sha*u_OQ*cos(zen)) - 2
		# Intermediate quantity u-path, p3, REF 2, Appendix A (A2), p. 656
		p3 = (c_isa**2*u_OQ**2*cos(zen)**2 + 2*c_isa*u_OQ*cos(zen) + 2)*exp(-c_isa*u_OQ*cos(zen)) - 2
		# Intermediate quantity u-path, p2, REF 2, Appendix A (A2), p. 656
		p2 = a_sha**-1*cos(zen*(1 - exp(-a_sha*u_OQ*cos(zen)) + ((16*p4*tan(zen)**2)/(9*pi*2*a_sha*R_T))))**-1
		# Intermediate quantity u-path, p1, REF 2, Appendix A (A2), p. 656
		p1 = c_isa**-1*cos(zen*(1 - exp(-c_isa*u_OQ*cos(zen)) + ((16*p3*tan(zen)**2)/(9*pi*2*c_isa*R_T))))**-1

		# Extinction of light along u-path from scatter at Q to observation at O, ksi1, REF 2, Appendix A (A2), p. 656
		ksi1 = exp(-N_m0*sig_m*(p1 + 11.778*K_am*p2)) 

		# Intermediate quantity s-path, f4, REF 2, Appendix A (A2), p. 657
		f4 = (a_sha**2*s_CQ**2*cos(psi)**2 + 2*a_sha*s_CQ*cos(psi) + 2)*exp(-a_sha*s_CQ*cos(psi)) - 2
		# Intermediate quantity s-path, f3, REF 2, Appendix A (A2), p. 657
		f3 = (c_isa**2*s_CQ**2*cos(psi)**2 + 2*c_isa*s_CQ*cos(psi) + 2)*exp(-c_isa*s_CQ*cos(psi)) - 2
		# Intermediate quantity s-path, f2, REF 2, Appendix A (A2), p. 657
		f2 = a_sha**-1*cos(psi*(1 - exp(-a_sha*s_CQ*cos(psi)) + ((16*f4*tan(psi)**2)/(9*pi*2*a_sha*R_T))))**-1
		# Intermediate quantity s-path, f1, REF 2, Appendix A (A2), p. 657
		f1 = c_isa**-1*cos(psi*(1 - exp(-c_isa*s_CQ*cos(psi)) + ((16*f3*tan(psi)**2)/(9*pi*2*c_isa*R_T))))**-1

		# Extinction of light along s-path from emission at C to scatter at Q, ksi2, REF 2, Appendix A (A2), p. 656
		ksi2 = exp(-N_m0*sig_m*(f1 + 11.778*K_am*f2))

		# Normalized emission function, I(psi), MODIFIED FROM REF 1, p. 13 (leaving out natural sky brightness)
		I_ne = 1/(2*pi)*(W_a*2*cos(psi) + W_b*0.554*psi**4 + W_c*sin(3*psi))

		# Illuminance per unit flux, i(psi,s), REF 2, Eq. 6, p. 644
		i_ps = I_ne*ksi2/s_CQ**2

		# Number density of gaseous component of atmosphere as function of altitude, N_m(h), REF 2, Eq. 10, p. 645
		N_m = N_m0*exp(-c_isa*h_Q)

		# Total integrated scattering cross-section, N_a*sig_a, REF 2, Eq. 12, p. 645 **REARRANGED**
		Na_x_siga = K_am*N_m*sig_m*11.11

		# Angular scattering function for molecular Rayleigh scattering, f_m, REF 2, Eq. 13, p. 646
		f_m = 3*(1 + cos(omega)**2)/(16*pi)

		# Angular scattering function for aerosol Mie scattering, f_a, REF 2, Eq. 14, p. 646
		if omega_deg >= 0.0 and omega_deg <= 10.0:
			f_a = 7.5*exp(-0.1249*omega_deg**2/(1 + 0.04996*omega_deg**2))
		elif omega_deg > 10.0 and omega_deg <= 124.0:
			f_a = 1.88*exp(-0.07226*omega_deg + 0.0002406*omega_deg**2)
		elif omega_deg > 124.0 and omega_deg <= 180.0:
			f_a = 0.025 + 0.015*sin((2.25*omega_deg - 369.0)*(pi/180))

		# Luminous flux per unit solid angle per unit upward flux (directly from source), S_d, REF 2, Eq. 5, p. 644
		S_d = (N_m*sig_m*f_m + Na_x_siga*f_a)*i_ps

		# Double scattering correction factor, D_S, REF 2, Eq. 20, p. 647
		D_S = 1 + N_m0*sig_m*(11.11*K_am*f2 + (f1/3))

		# Total illumance as a function of u, S(u), REF 2, Eq. 8, p. 645
		S_u = S_d*D_S
		
		# Integrand of propagation function, REF 2, Eq. 3, p. 644
		df_prop = S_u*ksi1*du
		total_sum += df_prop
		u_OQ += du
		lc += 1

	return total_sum

def array_to_geotiff(arr,outfilename,refVIIRS="20140901_20140930_75N180W_C.tif"):
	imdata = gdal.Open(refVIIRS)

	# Save out to a GeoTIFF
	[col,row] = arr.shape
	trans = imdata.GetGeoTransform()
	proj = imdata.GetProjection()
	nodatav = 0
	outfile = outfilename

	# Create the georeferenced file using the information from the VIIRS image
	outdriver = gdal.GetDriverByName("GTiff")
	outdata = outdriver.Create(str(outfile), row, col, 1, gdal.GDT_Float32)

	# Write data to the file, which is the kernel array in this example
	outdata.GetRasterBand(1).WriteArray(arr)

	# Set a no data value
	outdata.GetRasterBand(1).SetNoDataValue(nodatav)

	# Georeference the image
	outdata.SetGeoTransform(trans)

	# Write projection information
	outdata.SetProjection(proj)


def varrprint(varrval, varrtext):
	AddMessage("******************** {} :".format(varrtext))
	AddMessage("{} {}".format(varrval.shape,varrval.dtype))
	AddMessage("{} maximum: {}".format(varrtext, ma.maximum(varrval[~isnan(varrval)])))
	AddMessage("{} minimum: {}".format(varrtext, ma.minimum(varrval[~isnan(varrval)])))