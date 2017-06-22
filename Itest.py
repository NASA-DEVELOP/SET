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
import threading
# scipy has simpler implementation of Fast Fourier Transform (fft), ndimage.convolve is obsolete for our purposes with numpy.fft except as a comparison check
# https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.fftconvolve.html #scipy.signal.fftconvolve
from scipy import ndimage
import os.path
from matplotlib import pyplot as plt
import matplotlib.colors as colors
from matplotlib_scalebar.scalebar import ScaleBar
from osgeo import gdal

def main():
	# Print flag
	pflag = "verbose"
	kerneltiffpath = "kernelubreak30_highlat2.tif"
	if os.path.isfile(kerneltiffpath)==False:
		# Estimate the 2d propagation function and calc time and accuracy
		#bottom bottom_lat = 40.8797
		#top lat= 46.755666
		propagation_array1, time_1 = fsum_2d(pflag, 40.8797, 30.0)
		proparray_to_tiff(kerneltiffpath, propagation_array1)
		varrprint(propagation_array1,'propagation_array1', pflag)

		# propagation_array2, time_2 = fsum_2d(pflag, 40.8797, 10.0)
		# print "Time Factor Improvement!: {}".format(time_1/time_2)
		# varrprint(propagation_array2,'propagation_array2', pflag)
		# differencearray_perc = amax((abs(propagation_array1 - propagation_array2))/propagation_array1)
		# print "Accuracy Loss Factor!: {}".format(differencearray_perc)
		# varrprint(differencearray_perc, 'difference array perc', pflag)
		print "time for prop function ubreak 30"
		print time_1
	
	kerneldata = gdal.Open(kerneltiffpath)
	propagation_array1 = kerneldata.ReadAsArray()
	varrprint(propagation_array1, 'kernel',pflag)

	filein = "20140901_20140930_75N180W_C.tif"
	viirsraster = gdal.Open(filein)
	imagearr = viirsraster.ReadAsArray()
	varrprint(imagearr, 'VIIRS', pflag)

	######### Convert to float 32 for Fourier, scale, and round
	# falchi assumed natural sky brightness to be 174 micro cd/m^2 = 2.547e-11 watt/cm^2/steradian (at 555nm)
	# Not sure if this is correct scaling factor, I assume that this makes the output prop image in units of cd/m^2
	viirs_scaling_factor = 10**9
	imagearr *= viirs_scaling_factor
	propagation_array1 = float32(nan_to_num(propagation_array1))
	# generalized padding of kernel so that fft can run
	pad_left = (imagearr.shape[0] - propagation_array1.shape[0])//2
	pad_right = (imagearr.shape[0] - propagation_array1.shape[0])//2 + 1
	pad_up = (imagearr.shape[1] - propagation_array1.shape[1])//2
	pad_down = (imagearr.shape[1] - propagation_array1.shape[1])//2
	padded_prop = pad(propagation_array1,((pad_left,pad_right),(pad_up,pad_down)), 'constant', constant_values = 0)
	varrprint(padded_prop, 'kernel scaled and padded',pflag)
	varrprint(imagearr, 'VIIRSscaled', pflag)

	# subset = 50
	# prows = padded_prop.shape[0]
	# pcols = padded_prop.shape[1]
	# irows = padded_prop.shape[0]
	# icols = padded_prop.shape[1]

	# padded_prop = padded_prop[(prows//2)-subset:(prows//2)+subset, (pcols//2)-subset:(pcols//2)+subset]
	# imagearr = imagearr[(irows//2)-subset:(irows//2)+subset, (icols//2)-subset:(icols//2)+subset]
	# varrprint(padded_prop, "subsetted prop", pflag)
	# varrprint(imagearr, "subsetted viirs", pflag)

######################## Fourier Transform Method ################################
	def compare_arr(arr1, arr2, title1, title2, pixsize, norm1=True, norm2=True):
		scalebar = ScaleBar(pixsize)
		if norm1:
			plt.subplot(121),plt.imshow(arr1, norm = colors.LogNorm(), cmap = 'gray')
			plt.title(title1), plt.xticks([]), plt.yticks([])
			plt.gca().add_artist(scalebar)
		else:
			plt.subplot(121),plt.imshow(arr1, cmap = 'gray')
			plt.title(title1), plt.xticks([]), plt.yticks([])
			plt.gca().add_artist(scalebar)
		if norm2:
			plt.subplot(122),plt.imshow(arr2, norm = colors.LogNorm(), cmap = 'gray')
			plt.title(title2), plt.xticks([]), plt.yticks([])
		else:
			plt.subplot(122),plt.imshow(arr2, cmap = 'gray')
			plt.title(title2), plt.xticks([]), plt.yticks([])
		plt.show()

	np_dft_prop_im = fft.fft2(padded_prop)
	np_dft_kernel_shift = fft.fftshift(np_dft_prop_im)
	np_magnitude_spectrum = 20*log(abs(np_dft_kernel_shift))
	compare_arr(padded_prop, imagearr,'Relative Weights of Light Propogation (Convolution Kernel)', 'VIIRS Image', 462.7) #meters 462.7
	compare_arr(padded_prop, np_magnitude_spectrum,'Kernel', 'Fast Fourier Transformed Kernel',462.7, True, False) 
	varrprint(np_dft_prop_im, 'np_dft_kernel', pflag)

	np_dft_viirs_im = fft.fft2(imagearr)
	np_dft_viirs_shift = fft.fftshift(np_dft_viirs_im)
	np_magnitude_spectrum_viirs = 20*log(abs(np_dft_viirs_shift))
	compare_arr(imagearr, np_magnitude_spectrum_viirs,'VIIRS Image', 'Fast Fourier Transformed VIIRS', 462.7)
	varrprint(np_dft_viirs_im, 'np_dft_viirs', pflag)

	kernel_inv_shift = fft.ifftshift(np_dft_kernel_shift)
	viirs_inv_shift = fft.ifftshift(np_dft_viirs_shift)

	FFT_product_inverse = abs(fft.fftshift(fft.ifft2(kernel_inv_shift * viirs_inv_shift)))
	compare_arr(imagearr, FFT_product_inverse,'VIIRS Image', 'Product of FFT VIIRS and FFT Kernel: Artificial Light Propogation at Zenith', 462.7)
	varrprint(FFT_product_inverse, 'FFT_product_inverse', pflag)

	######## Comparison with Slow Convolution (Make sure to subset first) these give slightly different answers
	# apply kernel: https://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.ndimage.filters.convolve.html
	# filtered = ndimage.convolve(imagearr, padded_prop, mode='constant', cval = 0.0)
	# plt.subplot(121),plt.imshow(filtered, norm = colors.LogNorm(), cmap = 'gray')
	# plt.title('Slow Convolution Product'), plt.xticks([]), plt.yticks([])
	# plt.subplot(122),plt.imshow(FFT_product_inverse, norm = colors.LogNorm(), cmap = 'gray')
	# plt.title('Fast FFT Product'), plt.xticks([]), plt.yticks([])
	# plt.show()
	###############################################################################	
	
	proparray_to_geotiff(FFT_product_inverse)

# Function that creates 2d propagation function
def fsum_2d(pflag = 'verbose', latitude = 40.8797, ubr_arg = 10.0, zen_arg = 0.0, beta_arg = 0.0):
	# Input Variables
	print '**INPUTS**'

	#arbitrary radius and lat for testing purposes. Instead of R_teton to determine pixel should we use an array of radius of curvature?
	#bottom bottom_lat = 40.8797
	#top lat= 46.755666
	cent_lat_deg = latitude
	cent_lat = cent_lat_deg*pi/180
	p_deg = .0041666667
	p_rad = p_deg*pi/180

	# z, Zenith angle site, REF 2, Fig. 6, p. 648
	zen = zen_arg
	print 'z, Site zenith (deg): {}'.format(zen)

	# beta, Azimuth angle from line of sight to scatter from site, REF 2, Fig. 6, p. 648
	beta = beta_arg
	print 'beta, Relative azimuth line-of-sight to scatter (deg): {}'.format(beta)

	# ubr, Length of u for relaxing integration increment
	ubr = ubr_arg
	print 'ubr, Length of u for relaxing integration increment (km): {}'.format(ubr)

	# Earth radius at latitude 43.7904 (grand teton)
	R_teton = 6367.941 # km  (c)
	print 'R_teton, Radius of curvature of the Earth at Teton (km): {}'.format(R_teton)

	# Gaussian Earth radius of curvature, REF Wikipedia (https://en.wikipedia.org/wiki/Earth_radius)
	R_T = gauss_earth_curvature_radius(cent_lat)
	print 'R_T, Radius of curvature of the Earth (km): {}'.format(R_T)

	# Create latitude/longitude arrays
	source_lat, rltv_long, cent_row, cent_col = create_latlon_arrays(R_T, cent_lat, p_rad, pflag)

	# Distance from source (C) to observation site (O) along ellipsoid surface, REF 2, Fig. 6, p. 648
	# using haversine formula
	D_OC = 2*R_T*arcsin(sqrt(sin((source_lat - cent_lat)/2)**2 + cos(cent_lat)*cos(source_lat)*sin(rltv_long/2)**2))

	# Assignment of NaNs or null values outside of 200 km
	D_OC[D_OC > 201] = numpy.NaN

	# Check of D_OC shape after assigning NaNs outside of 200 km
	print "kernel heigth in pixels, trimmed: {}".format(D_OC.shape[0])
	print "kernel width in pixels, trimmed: {}".format(D_OC.shape[1])
	widthcenter = (D_OC.shape[1] + 1)//2
	heigthcenter = (D_OC.shape[0] +1)//2
	################################## reassignment of center value, need to use better method
	D_OC[cent_row,cent_col] = .01
	varrprint(D_OC,'D_OC', pflag)
	# print "center area ##################################################################"
	# print D_OC[heigthcenter-2:heigthcenter+1,widthcenter-2:widthcenter+1]
	# Earth angle from source to site, REF 3, p. 308\
	Chi = D_OC/R_T
	varrprint(Chi,'Chi', pflag)
	
	# u0, shortest scattering distance based on curvature of the Earth, REF 2, Eq. 21, p. 647
	u0 = 2*R_T*sin(Chi/2)**2/(sin(zen)*cos(beta)*sin(Chi)+cos(zen)*cos(Chi)) #km
	varrprint(u0,'u0', pflag)
	
	# l, Direct line of sight distance between source and observations site, REF 2, Appendix A (A1), p. 656
	# l_OC and D_OC are similar as expected
	l_OC = sqrt(4*R_T**2*sin(Chi/2)**2) # km

	print l_OC[heigthcenter-2:heigthcenter+1,widthcenter-2:widthcenter+1]
	# q1, Intermediate quantity, REF 2, Appendix A (A1), p. 656, **WITH CORRECTION FROM REF 3, eq. 6, p. 308**
	q1 = R_T*(sin(Chi)*sin(zen)*cos(beta) + cos(Chi)*cos(zen) - cos(zen)) # km

	print q1[heigthcenter-2:heigthcenter+1,widthcenter-2:widthcenter+1]
	# theta, elevation angle of scatter above source from site (QOC), REF 2, Appendix A (A1), p. 656
	theta = arccos(q1/l_OC) # radians
	
	print theta[heigthcenter-2:heigthcenter+1,widthcenter-2:widthcenter+1]
	# Get left arrays to cut processing time in half
	Chileft = Chi[0:,0:widthcenter]
	u0left = u0[0:,0:widthcenter]
	l_OCleft = l_OC[0:,0:widthcenter]
	thetaleft = theta[0:,0:widthcenter]

	#test array subsets to reduce processing time
	# Chileft = Chi[427:432,530:widthcenter]
	# u0left = u0[427:432,530:widthcenter]
	# l_OCleft = l_OC[427:432,530:widthcenter]
	# thetaleft = theta[427:432,530:widthcenter]

	#container for Propogation array
	PropSumArrayleft = zeros_like(l_OCleft)

	# ################## Threading # currently does no better than no threading
	# # print "Time in seconds for 400 iterations, with threads"

	# # start = time.time()

	# # threads = []
	# # def helper(p, p2p, R_T, c, u, l, t):
	# # 	p[...] = fsum(p2p, R_T, c, u, l, t)

	# # for p, c, u, l, t in itertools.izip(nditer(PropSumArray, op_flags=['readwrite']),nditer(Chi, op_flags=['readwrite']),nditer(u0, op_flags=['readwrite']), nditer(l_OC, op_flags=['readwrite']),nditer(theta, op_flags=['readwrite'])):
	# #     threads.append(threading.Thread(target=helper(p, p2p, R_T, c, u, l ,t), args=(p, p2p, R_T, c, u, l, t)))

	# # count = 0
	# # for thread in threads:
	# # 	if (count % 30 == 0 and count != 0):
	# # 		time.sleep(.6) #time for previous threads to finish so the next round doesn't interfere with the previous
	# # 	thread.start()
	# # 	count += 1

	# # end = time.time()
	# # print (end-start)

	# ############################
	
	# ################### NoThreading #####################
	print "Time for iterations, no threads"
	start = time.time()

	# 2d iteration for integrating from u0 to infinity to create propagation function for each element
	for p,c,u,l,t in itertools.izip(nditer(PropSumArrayleft, op_flags=['readwrite']),nditer(Chileft, op_flags=['readwrite']),nditer(u0left, op_flags=['readwrite']), nditer(l_OCleft, op_flags=['readwrite']),nditer(thetaleft, op_flags=['readwrite'])):
		p[...] = fsum_single(R_T, c, u, l, t, zen, beta, ubr)
	end = time.time()
	time_sec = end-start
	PropSumArrayright = fliplr(PropSumArrayleft[:,1:])

	# Complete 2d propagation function
	PropSumArray = hstack((PropSumArrayleft, PropSumArrayright))
	print "prop sum array hstacked"
	print PropSumArray
	return PropSumArray, time_sec

# Function to calculate Gaussian Earth radius of curvature as a function of latitude
def gauss_earth_curvature_radius(center_lat):
	# Earth ellipse semi-major orbit, REF Wikipedia (https://en.wikipedia.org/wiki/Earth_radius)
	R_equator = 6378.1370 # km (a)

	# Earth ellipse semi-minor axis, REF Wikipedia (https://en.wikipedia.org/wiki/Earth_radius)
	R_polar = 6356.7523142 # km (b)

	# Gaussian Earth radius of curvature, REF Wikipedia (https://en.wikipedia.org/wiki/Earth_radius)
	R_curve = ((R_equator**2)*R_polar)/((R_equator*cos(center_lat))**2 + (R_polar*sin(center_lat))**2)
	
	# NOT USED PRESENTLY: Directional Earth radius of curvature, REF Wikipedia (https://en.wikipedia.org/wiki/Earth_radius)
	# R_T = (R_polar*R_equator**2)/((R_equator*cos((cent_lat_rad+rel_lat_rad/2)))**2+(R_polar*sin((cent_lat_rad+rel_lat_rad/2)))**2)

	return R_curve

# Function does initial array sizing and create latitude/longitude arrays
def create_latlon_arrays(R_curve, center_lat, pix_rad, pf):
	p_h = R_curve*pix_rad
	p_w = cos(center_lat)*R_curve*pix_rad

	# Calculate dimensions for kernel
	kernel_cols = int(round(400/p_w))
	kernel_rows = int(round(400/p_h))
	if kernel_cols%2 == 0:
		kernel_cols += 1 
	if kernel_rows%2 == 0:
		kernel_rows += 1
	print "kernel width in columns, untrimmed: {}".format(kernel_cols)
	print "kernel height in rows, untrimmed: {}".format(kernel_rows)

	# Create vectors of column and row counts
	col_count = array(range(kernel_cols))
	row_count = array(range(kernel_rows))

	# Find center column, center row
	center_col = int((kernel_cols - 1)/2)
	center_row = int((kernel_rows - 1)/2)
	print center_col
	print center_row

	# Create vecotrs of relative longitude and latitudes (in radians)
	rel_long_vec = array(pix_rad*(col_count - center_col))
	print rel_long_vec.shape
	rel_lat_vec = -pix_rad*(row_count - center_row)
	print rel_lat_vec.shape
	print rel_long_vec
	print rel_lat_vec
	print rel_long_vec[center_col]
	print rel_lat_vec[center_row] 

	rel_long = tile(rel_long_vec,(kernel_rows,1))
	rel_lat = transpose(tile(rel_lat_vec,(kernel_cols,1)))
	print rel_long
	varrprint(rel_long,'rel_long',pf)
	print rel_lat
	varrprint(rel_lat,'rel_lat',pf)

	src_lat = rel_lat + center_lat
	print src_lat
	varrprint(src_lat,'src_lat',pf)

	return src_lat, rel_long, center_row, center_col

# Function that takes elements of the arrays of D_OC, Chi, etc. as arguemnts. 
def fsum_single(R_T, Chi, u0, l_OC, theta, zen_farg, beta_farg, ubrk_farg, K_am_arg = 1.0, del_u_farg = .2):
	if isnan(l_OC):
		return nan

	# Scattering distance increment
	del_u0 = del_u_farg
	# print 'delta_u, Scattering distance increment for finite integration over u (km): {}'.format(del_u)

	# zenith angle at observation site
	zen = zen_farg # radians

	# horizontal angle at observation site from line-of-sight (OC) to line of incoming scattered light (OQ)
	beta = beta_farg #radians

	# length of u, scattering distance, at which the integration increment is relaxed
	ubreak = ubrk_farg #km

	# K - Parameter for relative importance of aerosols to molecules, REF 1, p. 10
	K_am = K_am_arg # unitless ratio

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

	# loop counter
	lc = 1
	
	while u_OQ < 30.0: # df_prop > stability_limit*total_sum:
		## START OF "u" LOOP

		if u_OQ < ubreak:
			del_u = del_u0
		else:
			del_u = del_u0*10.0

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
		lc += 1
		# print lc,

	return total_sum

def proparray_to_geotiff(convolved_array, referenceVIIRS="C:/MonthlyViirs2015/20140901_20140930_75N180W_C.tif", outfilename= 'C:/outputkerneltiffs/referenced.tif'):
	imdata = gdal.Open(referenceVIIRS)

	# Save out to a GeoTiff
	arr = convolved_array
	# First of all, gather some information from VIIRS image
	[cols,rows] = arr.shape
	trans       = imdata.GetGeoTransform()
	proj        = imdata.GetProjection()
	nodatav     = 0
	outfile     = outfilename

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

if __name__ == "__main__":
	main()