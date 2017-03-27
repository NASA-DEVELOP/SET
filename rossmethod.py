# finsih distance function, figure out degree radians conversion, which variable needs which units.

# REFERENCES
# (1) Falchi, F., P. Cinzano, D. Duriscoe, C.C.M. Kyba, C.D. Elvidge, K. Baugh, B.A. Portnov, N.A. Rybnikova
#       and R. Furgoni, 2016. The new workd atlas of artificial night sky brightness. Sci. Adv. 2.
# (2) Cinzano, P., F. Falchi, C.D. Elvidge and  K.E. Baugh, 2000. The artificial night sky brightness mapped
#       from DMSP satellite Operational Linescan System measurements. Mon. Not. R. Astron. Soc. 318.
# (3) Garstang, R.H., 1989. Night-sky brightness at observatories and sites. Pub. Astron. Soc. Pac. 101.

from __future__ import division
import numpy
from numpy import *
import itertools
import time
import archook
archook.get_arcpy()
import arcpy
import threading

class BrightP2P(object):
	def __init__(self, K_am):
		# K - Parameter for relative importance of aerosols to molecules, REF 1, p. 10
		self.K_am = K_am # unitless ratio

	def set_zenith(self, z):
		self.z = z

	def set_beta(self, beta):
		self.beta = beta

	def set_del_u(self, del_u):
		self.del_u = del_u

	def set_radius_curvature(self):
		self.R_T = (self.R_polar*self.R_equator**2)/((self.R_equator*cos((self.lat_O+self.lat_C/2)))**2+(self.R_polar*sin((self.lat_O+self.lat_C/2)))**2)

	def get_z_deg(self):
		return self.z*180/pi

	def get_beta_deg(self):
		return self.beta*180/pi

## Functions ##

# funciton to slice arrays down the middle, then stitch back together after adding pixel increment element wise
def slice_n_stack(arr, pix_measure, central_index):
	arr_left = arr.T[0:central_index]
	
	arr_right = arr.T[central_index+1:]
	
	p_incr = pix_measure
	for column in arr_left:
		column -= p_incr
		p_incr += pix_measure
	arr_left = arr_left[ : : -1]
	print "shape left: {}".format(arr_left.shape)
	p_incr = pix_measure
	for column in arr_right:
		column += p_incr
		p_incr += pix_measure
	print "shape right: {}".format(arr_right.shape)
	centercolumn = arr.T[central_index:central_index+1]
	
	result = hstack((arr_left.T, centercolumn.T, arr_right.T))
	return result

# Function that takes elements of the arrays of D_Oc, Chi , etc. as arguemnts. 
# Arguments named the same as arrays for laziness

def fsum(p2p, R_T, Chi, u0, l_OC, theta):
	
	if isnan(l_OC):
		return nan
	# Unpack point-to-point variables. Class may not be necessary anymore
	del_u = p2p.del_u
	K_am = p2p.K_am

	#*Constants*
	# N_m,0 - Molecular density at sea level, REF 2, p. 645
	N_m0 = 2.55e19 # cm^-3

	# c - Inverse scale altitude, REF. 2, p. 645
	c_isa = 0.104 # km^-1

	# sigma_m - Integrated Rayleigh scattering cross-section in the V-band, REF 2, p. 646
	sig_m = 1.136e-26 # cm^-2*sr^-1

	# Falchi best fit parameters for normalzied emssion
	# Weights for angular distributions. REF 1, p. 21
	W_a = 1.9e-3
	W_b = 5.2e-4
	W_c = 7.6e-5
	# light loss per night hour after midnight (negative hours before midnight)
	d_ll = -0.045
	#*End Constants*
	
	#containers for variables that update in loop
	u_OQ = u0 # km
	total_sum = 0
	df_prop = 1

	#Total Propogation stable to 3 significant figures
	stability_limit = 0.001
	
	while df_prop > stability_limit*total_sum:
		## START OF "u" LOOP
		# s, Distance from source to scattering CQ, REF 2, Appendix A (A1), p. 656
		# equation is wrong in Ref 2 (Cinzano). Changed angle from Chi to theta according to Ref 3, p. 308, Equation 7 (Garstang)
		s_CQ = sqrt((u_OQ - l_OC)**2 + 4*u_OQ*l_OC*sin(theta/2)**2) # km
		
		# h, Height of scattering above Earth reference surface Q, REF 2, Appendix A (A1), p. 656
		h_Q = R_T*(sqrt(1 + (u_OQ**2 + 2*u_OQ*R_T*cos(zen))/R_T**2) - 1) # km
		

		# phi, Elevation angle of emission from C to Q over line of sight from C to O (QCO), REF 2, Appendix A (A1), p. 656
		phi = arccos((l_OC**2 + s_CQ**2 - u_OQ**2)/(2*l_OC*s_CQ)) # radians
		phi_deg = phi*180/pi # degrees
		

		# q3, Intermediate quantity, REF 2, Appendix A (A1), p. 656 
		q3 = u_OQ*cos(zen)*cos(Chi) - 2*R_T*sin(Chi/2)**2 # km
		

		# q2, Intermediate quantity (q2=q3 if zenith angle, z=0), REF 2, Appendix A (A1), p. 656
		q2 = u_OQ*sin(zen)*cos(beta)*sin(Chi) + q3 # km
		

		# Psi, emission angle from source, REF 2, Appendix A (A1), p. 656
		Psi = arccos(q2/s_CQ) # radians
		Psi_deg = Psi*180/pi # degrees
		

		# omega, scattering angle at Q, REF 2, Appendix A (A1), p. 656
		omega = theta + phi # radians
		omega_deg = omega*180/pi # degrees
		

		
		# a, Scale height of aerosols, REF 2, p. 646
		a_sha = 0.657 + 0.059*K_am # km^-1
		

		# p4, Intermediate quantity u-path, REF 2, Appendix A (A2), p. 656
		p4 = (a_sha**2*u_OQ**2*cos(zen)**2 + 2*a_sha*u_OQ*cos(zen) + 2)*exp(-a_sha*u_OQ*cos(zen)) - 2
		

		# p3, Intermediate quantity u-path, REF 2, Appendix A (A2), p. 656
		p3 = (c_isa**2*u_OQ**2*cos(zen)**2 + 2*c_isa*u_OQ*cos(zen) + 2)*exp(-c_isa*u_OQ*cos(zen)) - 2
		

		# p2, Intermediate quantity u-path, REF 2, Appendix A (A2), p. 656
		p2 = a_sha**-1*cos(zen*(1 - exp(-a_sha*u_OQ*cos(zen)) + ((16*p4*tan(zen)**2)/(9*pi*2*a_sha*R_T))))**-1
		

		# p1, Intermediate quantity u-path, REF 2, Appendix A (A2), p. 656
		p1 = c_isa**-1*cos(zen*(1 - exp(-c_isa*u_OQ*cos(zen)) + ((16*p3*tan(zen)**2)/(9*pi*2*c_isa*R_T))))**-1
		

		# ksi1, Extinction of light along u-path from scatter at Q to observation at O, REF 2, Appendix A (A2), p. 656
		ksi1 = exp(-N_m0*sig_m*(p1 + 11.778*K_am*p2)) 
		

		# f4, Intermediate quantity s-path, REF 2, Appendix A (A2), p. 657
		f4 = (a_sha**2*s_CQ**2*cos(Psi)**2 + 2*a_sha*s_CQ*cos(Psi) + 2)*exp(-a_sha*s_CQ*cos(Psi)) - 2
		

		# f3, Intermediate quantity s-path, REF 2, Appendix A (A2), p. 657
		f3 = (c_isa**2*s_CQ**2*cos(Psi)**2 + 2*c_isa*s_CQ*cos(Psi) + 2)*exp(-c_isa*s_CQ*cos(Psi)) - 2
		

		# f2, Intermediate quantity s-path, REF 2, Appendix A (A2), p. 657
		f2 = a_sha**-1*cos(Psi*(1 - exp(-a_sha*s_CQ*cos(Psi)) + ((16*f4*tan(Psi)**2)/(9*pi*2*a_sha*R_T))))**-1
		

		# f1, Intermediate quantity s-path, REF 2, Appendix A (A2), p. 657
		f1 = c_isa**-1*cos(Psi*(1 - exp(-c_isa*s_CQ*cos(Psi)) + ((16*f3*tan(Psi)**2)/(9*pi*2*c_isa*R_T))))**-1
		

		# ksi2, Extinction of light along s-path from emission at C to scatter at Q, REF 2, Appendix A (A2), p. 656
		ksi2 = exp(-N_m0*sig_m*(f1 + 11.778*K_am*f2)) 
		

		# I(Psi), Normalized emission function, MODIFIED FROM REF 1, p. 13 (leaving out natural sky brightness)
		I_ne = 1/(2*pi)*(W_a*2*cos(Psi) + W_b*0.554*Psi**4 + W_c*sin(3*Psi))
		

		# i(Psi,s), Illuminance per unit flux, REF 2, Eq. 6, p. 644
		i_ps = I_ne*ksi2/s_CQ**2
		

		# N_m(h), Number density of gaseous component of atmosphere as function of altitude, h, REF 2, Eq. 10, p. 645
		N_m = N_m0*exp(-c_isa*h_Q)
		

		# N_a*sigma_a, Total integrated scattering cross-section, REF 2, Eq. 12, p. 645 **REARRANGED**
		Na_x_siga = K_am*N_m*sig_m*11.11
		

		# f_m, Angular scattering function for molecular Rayleigh scattering, REF 2, Eq. 13, p. 646
		f_m = 3*(1 + cos(omega)**2)/(16*pi)

		# f_a, Angular scattering function for aerosol Mie scattering, REF 2, Eq. 14, p. 646
		if omega_deg >= 0.0 and omega_deg <= 10.0:
			f_a = 7.5*exp(-0.1249*omega_deg**2/(1 + 0.04996*omega_deg**2))
		elif omega_deg > 10.0 and omega_deg <= 124.0:
			f_a = 1.88*exp(-0.07226*omega_deg + 0.0002406*omega_deg**2)
		elif omega_deg > 124.0 and omega_deg <= 180.0:
			f_a = 0.025 + 0.015*sin((2.25*omega_deg - 369.0)*(pi/180))
		

		# S_d, Luminous flux per unit solid angle per unit upward flux (directly from source), REF 2, Eq. 5, p. 644
		S_d = (N_m*sig_m*f_m + Na_x_siga*f_a)*i_ps
		

		# D_S, Double scattering correction factor, REF 2, Eq. 20, p. 647
		D_S = 1 + N_m0*sig_m*(11.11*K_am*f2 + (f1/3))
		

		# S_u, Total illumance as a function of u, REF 2, Eq. 8, p. 645
		S_u = S_d*D_S
		
		df_prop = S_u*ksi1*del_u
		# integrand of propogation function, REF 2, Eq. 3, p. 644
		total_sum = df_prop + total_sum
		u_OQ += del_u

	return total_sum

# array variable check print function
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

# Input Variables
print '**INPUTS**'

# Print flag
pflag = 'verbose'

# Initialize BrightP2P instance that has point-to-point properties
p2p = BrightP2P(1.0)

# z, Zenith angle site, REF 2, Fig. 6, p. 648
p2p.set_zenith(0.0) # deg
print 'z, Site zenith (deg): {}'.format(p2p.get_z_deg())

# beta, Azimuth angle from line of sight to scatter from site, REF 2, Fig. 6, p. 648
p2p.set_beta(0.0) #deg
print 'beta, Relative azimuth line-of-sight to scatter (deg): {}'.format(p2p.get_beta_deg())

# Scattering distance increment
p2p.set_del_u(0.2) #km
print 'delta_u, Scattering distance increment for finite integration over u (km): {}'.format(p2p.del_u)

# # Earth radius of curvature, REF Wikipedia (https://en.wikipedia.org/wiki/Earth_radius)
# p2p.set_radius_curvature() #km
# print 'R_T, Radius of curvature of the Earth (km): {}'.format(p2p.R_T)

################################################################# Moving Window Kernel
#single case array build
#import arcpy

# Earth ellipse semi-major orbit, REF Wikipedia (https://en.wikipedia.org/wiki/Earth_radius)
R_equator = 6378.1370 # km (a)

# Earth ellipse semi-minor axis, REF Wikipedia (https://en.wikipedia.org/wiki/Earth_radius)
R_polar = 6356.7523142 # km (b)

# Earth radius at latitude 43.7904 (grand teton)
R_teton = 6367.941 # km  (c)

zen = p2p.z

beta = p2p.beta

#arbitrary radius and lat for testing purposes. Instead of R_teton to determine pixel should we use an array of radius of curvature?
#bottom cent_lat = 40.8797
#top
cent_lat = 46.755666
p_deg = .0041666667

p_h = R_teton*p_deg*pi/180
p_w = cos(cent_lat*pi/180)*R_teton*p_deg*pi/180
cent_lat_km = cent_lat*R_teton*pi/180
print "cent lat"
print cent_lat_km
#pixel units

# Calculate dimensions for kernel
kernel_cols = int(round(400/p_w))
kernel_rows = int(round(400/p_h))
if kernel_cols%2 == 0:
	kernel_cols += 1 
if kernel_rows%2 == 0:
	kernel_rows += 1
print "kernel width in columns, untrimmed: {}".format(kernel_cols)
print "kernel height in rows, untrimmed: {}".format(kernel_rows)
col_count = array(range(kernel_cols))
row_count = array(range(kernel_rows))
center_col = int((kernel_cols - 1)/2)
center_row = int((kernel_rows - 1)/2)
print center_col
print center_row
rel_long_vec = p_deg*(col_count - center_col)
print rel_long_vec.shape
rel_lat_vec = -p_deg*(row_count - center_row)
print rel_lat_vec.shape
print rel_long_vec
print rel_lat_vec
print rel_long_vec[center_col]
print rel_lat_vec[center_row] 

rel_long0 = tile(rel_long_vec,(kernel_rows,1))
rel_lat0 = transpose(tile(rel_lat_vec,(kernel_cols,1)))
print rel_long0
varrprint(rel_long0,'rel_long0',pflag)
#are these minimums and maximums correct for rel_lat0? they look to be the same as rel_long0. Rel_long0 looks shaped like latitude
print rel_lat0
varrprint(rel_lat0,'rel_lat0',pflag)

#should cent_lat be appended here?
lat_array = rel_lat0 + cent_lat
print lat_array
varrprint(lat_array,'lat_array',pflag)

# #made relative arrays the same size, easier, few extra distance calculations will be trimmed
# rel_long = zeros((int(kernel_rows), int(kernel_cols)))
# # print "shape long initial: {}".format(rel_long.shape)
# rel_lat = ones((int(kernel_rows), int(kernel_cols)))*cent_lat_km
# # print "shape lat initial: {}".format(rel_long.shape)
# centralcolumnindex_long = int(round(rel_long.shape[1]/2))-1
# centralrowindex_lat = int(round(rel_lat.shape[0]/2))-1
# print "Central column: {}".format(centralcolumnindex_long)
# print "Central row: {}".format(centralrowindex_lat)
# rel_long_km = slice_n_stack(rel_long, p_h, centralcolumnindex_long)
# # print "shape long km: {}".format(rel_long_km.shape)
# rel_lat_km = slice_n_stack(rel_lat.T, p_w, centralrowindex_lat)
# # print "shape lat km: {}".format(rel_long_km.shape)
# varrprint(rel_long_km,'relative longitude array', pflag)
# varrprint(rel_lat_km,'relative latitude array', pflag)

# rel_lat_km = rel_lat_km.T
# varrprint(rel_lat_km,'relative latitude array *T*', pflag)

# print"******************** Relative Longitude (KM) Array"
# print rel_long_km[0:1,centralcolumnindex_long:centralcolumnindex_long+20]


# print"******************** Relative Lat (KM) Array"
# print rel_lat_km[centralrowindex_lat:centralrowindex_lat+20,0:1]
# # convert to radians for haversine formula
rel_long_rad = rel_long_km/R_teton
rel_lat_rad = rel_lat_km/R_teton
cent_lat_rad = cent_lat*pi/180
# print "shape long: {}".format(rel_long_rad.shape)
# print "shape lat: {}".format(rel_lat_rad.shape)
# Distance from source (C) to observation site (O) along ellipsoid surface, REF 2, Fig. 6, p. 648
D_OC = 2*R_teton*arcsin(sqrt(sin((cent_lat_rad-rel_lat_rad)/2)**2 + cos(rel_lat_rad)*cos(cent_lat_rad)*sin((rel_long_rad-0)/2)**2))
# print 'D_OC array'
# print "Center Distance, is close too but not equal to zero"
# print D_OC[centralrowindex_lat:centralrowindex_lat+1, centralcolumnindex_long:centralcolumnindex_long+1]
# # Radius of curvature
# # R_T = (R_polar*R_equator**2)/((R_equator*cos((cent_lat_rad+rel_lat_rad/2)))**2+(R_polar*sin((cent_lat_rad+rel_lat_rad/2)))**2)
# R_T = 6367.941
# # Can't figure out how to remove nans from this array, will go with "good option for now"
# # for x,y in itertools.izip(nditer(D_OC, op_flags=['readwrite']), nditer(R_T, op_flags=['readwrite'])):
# # 	if x > 201:
# # 		y[...] = numpy.NaN

# D_OC[D_OC > 201] = numpy.NaN

# #functions to prune all columns with NaNs. Need to figure out elemen wiseMay also need to prune corresponding elements in lat and longitude arrays
# #D_OC = D_OC[:, ~isnan(D_OC).any(axis=0)]


# ##causes radius_kernel to be empty
# #radius_kernel = radius_kernel[:, ~isnan(radius_kernel).any(axis=0)]
# print "kernel width in pixels, trimmed: {}".format(D_OC.shape[0])

# print "kernel height in pixels, trimmed: {}".format(D_OC.shape[1])

# print"******************** D_OC Array"
# ##################################reassignment of center value, need to use better method
# D_OC[630:631,630:631] =.01
# varrprint(D_OC,'D_OC', pflag)

# # Earth angle from source to site, REF 3, p. 308
# Chi = D_OC/R_T
# varrprint(Chi,'Chi', pflag)

# #u0, shortest scattering distance based on curvature of the Earth, REF 2, Eq. 21, p. 647
# u0 = 2*R_T*sin(Chi/2)**2/(sin(zen)*cos(beta)*sin(Chi)+cos(zen)*cos(Chi)) #km
# varrprint(u0,'u0', pflag)

# # l, Direct line of sight distance between source and observations site, REF 2, Appendix A (A1), p. 656
# # L_OC and D_OC are similar as expected
# l_OC = sqrt(4*R_T**2*sin(Chi/2)**2) # km
# varrprint(l_OC,'l_OC', pflag)

# # q1, Intermediate quantity, REF 2, Appendix A (A1), p. 656, **WITH CORRECTION FROM REF 3, eq. 6, p. 308**
# q1 = R_T*(sin(Chi)*sin(zen)*cos(beta) + cos(Chi)*cos(zen) - cos(zen)) # km
# varrprint(q1,'q1', pflag)

# # theta, elevation angle of scatter above source from site (QOC), REF 2, Appendix A (A1), p. 656
# theta = arccos(q1/l_OC) # radians
# varrprint(theta,'theta', pflag)


# Chileft = Chi[0:,0:632]
# u0left = u0[0:,0:632]
# l_OCleft = l_OC[0:,0:632]
# thetaleft = theta[0:,0:632]

# # print"******************** D_OCleft Array Subsetted"
# # print D_OCleft
# # print"******************** Chileft Array Sub"
# # print Chileft
# # print"******************** u0left Array Sub"
# # print u0left
# # print"******************** l_OCleft Array Sub"
# # print l_OCleft
# # print"******************** thetaleft Array Sub"
# # print thetaleft
# # https://docs.scipy.org/doc/numpy/reference/arrays.nditer.html Iterator-Allocated output Arrays
# PropSumArrayleft = zeros_like(l_OCleft)

# # ################## Threading # currently does no better than no threading
# # # print "Time in seconds for 400 iterations, with threads"

# # # start = time.time()

# # # threads = []
# # # def helper(p, p2p, R_T, c, u, l, t):
# # # 	p[...] = fsum(p2p, R_T, c, u, l, t)

# # # for p, c, u, l, t in itertools.izip(nditer(PropSumArray, op_flags=['readwrite']),nditer(Chi, op_flags=['readwrite']),nditer(u0, op_flags=['readwrite']), nditer(l_OC, op_flags=['readwrite']),nditer(theta, op_flags=['readwrite'])):
# # #     threads.append(threading.Thread(target=helper(p, p2p, R_T, c, u, l ,t), args=(p, p2p, R_T, c, u, l, t)))

# # # count = 0
# # # for thread in threads:
# # # 	if (count % 30 == 0 and count != 0):
# # # 		time.sleep(.6) #time for previous threads to finish so the next round doesn't interfere with the previous
# # # 	thread.start()
# # # 	count += 1

# # # end = time.time()
# # # print (end-start)

# # ############################


# # ################### NoThreading
# print "Time in seconds for 400 iterations, no threads"
# start = time.time()
# # for p,c,u,l,t in itertools.izip(nditer(PropSumArrayleft, op_flags=['readwrite']),nditer(Chileft, op_flags=['readwrite']),nditer(u0left, op_flags=['readwrite']), nditer(l_OCleft, op_flags=['readwrite']),nditer(thetaleft, op_flags=['readwrite'])):
# # 	p[...] = fsum(p2p, R_T, c, u, l, t)


# end = time.time()
# print (end-start)
# # PropSumArrayright = fliplr(PropSumArrayleft[:,1:])
# # PropSumArray = hstack((PropSumArrayleft, PropSumArrayright))
# # D_OCright = fliplr(D_OCleft[:,1:])
# # D_OCarr = hstack((D_OCleft, D_OCright))
# # print "array shapes"
# # print PropSumArray
# # print D_OCarr
# long_deg = rel_long_rad*180/pi
# lat_deg = rel_lat_rad*180/pi

# ########################### Array to Raster
# # filein = "C:/VIIRS_processing/Clipped Rasters.gdb/VIIRS_2014_06"
# # myRaster = arcpy.Raster(filein)
# # arcpy.env.overwriteOutput = True
# # arcpy.env.outputCoordinateSystem = filein
# # arcpy.env.cellSize = filein

# # latitude = 44.384922
# # longitude = -110.650200

# # lower_left = arcpy.Point(longitude,latitude)

# # x_size = cos(cent_lat*pi/180)*p_deg
# # y_size = p_deg

# # Kernel = arcpy.NumPyArrayToRaster(PropSumArray, lower_left, x_size, y_size, nan)
# # output = "C:/ArtificialBrightness/BigKerneltest2.tif"
# # Kernel.save(output)
# #####################

# # print "*************************Propogation Array*******************************"
# savetxt("timetest.txt", PropSumArrayleft, fmt= "%.6e", delimiter= ',', newline=';')




#didn't end up using ellipsoidal calculation of distance, we used great circle. 
#if we decide we can use this guy's code for ellipsoid calc

# **) MIT License
#
# Copyright (C) 2016-2017 -- mrJean1@Gmail.com
#
# Permission is hereby granted, free of charge, to any person obtaining a
# copy of this software and associated documentation files (the "Software"),
# to deal in the Software without restriction, including without limitation
# the rights to use, copy, modify, merge, publish, distribute, sublicense,
# and/or sell copies of the Software, and to permit persons to whom the
# Software is furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included
# in all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
# OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL
# THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR
# OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
# ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
# OTHER DEALINGS IN THE SOFTWARE.