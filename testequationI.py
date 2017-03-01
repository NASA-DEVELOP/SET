# In the words of Cinzano 2000, here is what we are trying to calculate: 
# "The light pollution propagation function f depends in general on the geometrical disposition (the altitude of the site and the area, and their mutual distance), on the atmospheric
# distribution of molecules and aerosols and their optical characteristics in the chosen photometric band, on the shape of the emission function of the source, and on the direction 
# of the sky observed." 

# Our goal is to get a rough calculation of the illuminance per unit flux along the line of sight of an observer looking straight up (at the zenith), 
# using Dr. Ross's estimation of the light propogation functions, which are sketched out in Cinzano 2000 and Falchi 2016.

# We will first just work on translating these algebraic equations to python code for a given emission angle. We will be using a package called numpy, which has the capability to work with arrays. 
# Arrays are rasters that can't be loaded up into geoprocessing software. They are quick to process and are a simple data format for data that occurs in a matrix.

#Second we will work on integrating each of these algebraic equations to calculate the total amount of radiance occuring from the possible range of emission angles (psi).

#Third, we will apply the integral equations over all "pixels" in the array using a nested for loop (for loop within a for loop).

# Here are our equations to describe upward light emission

# I(Psi)_A = (1/(2*pi))*2*cos(Psi)   # This is the Lambertian circular distribution for stuff reflected from the ground
# I(Psi)_B = (1/(2*pi))*0.554*Psi^4  # This is the horizon lighting model poorly shielded street lighting
# I(Psi)_C = (1/(2*pi))*sin(3*Psi)   # This is the rose petal function that is max at 30 degrees. It is for the intermediate light

# REFERENCES
# (1) Falchi, F., P. Cinzano, D. Duriscoe, C.C.M. Kyba, C.D. Elvidge, K. Baugh, B.A. Portnov, N.A. Rybnikova
#       and R. Furgoni, 2016. The new workd atlas of artificial night sky brightness. Sci. Adv. 2.
# (2) Cinzano, P., F. Falchi, C.D. Elvidge and  K.E. Baugh, 2000. The artificial night sky brightness mapped
#       from DMSP satellite Operational Linescan System measurements. Mon. Not. R. Astron. Soc. 318.
# (3) Garstang, R.H., 1989. Night-sky brightness at observatories and sites. Pub. Astron. Soc. Pac. 101.

import numpy
from numpy import *

class BrightP2P(object):
	# Earth ellipse semi-major orbit, REF Wikipedia (https://en.wikipedia.org/wiki/Earth_radius)
	R_equator = 6378.1370 # km (a)

	# Earth ellipse semi-minor axis, REF Wikipedia (https://en.wikipedia.org/wiki/Earth_radius)
	R_polar = 6356.7523142 # km (b)

	def set_lat(self, lat0):
		self.lat0 = lat0

	def set_radius_curvature(self):
		self.R_T = (self.R_polar*self.R_equator**2)/((self.R_equator*cos(self.lat0))**2+(self.R_polar*sin(self.lat0))**2)

	def set_distance(self, D):
		self.D = D

	def set_zenith(self, z):
		self.z = z

	def set_beta(self, beta):
		self.beta = beta

	def set_del_u(self, del_u):
		self.del_u = del_u

	def get_lat_dd(self):
		return self.lat0*180/pi

	def get_z_deg(self):
		return self.z*180/pi

	def get_beta_deg(self):
		return self.beta*180/pi

# Input Variables
print '**INPUTS**'

# Initialize BrightP2P instance that has point-to-point properties
p2p = BrightP2P()

# Latitude of observer
p2p.set_lat(43.5*pi/180) # degrees
print 'Latitude (deg): {}'.format(p2p.get_lat_dd())

# D, Distance from source (C) to observation site (O) along ellipsoid surface, REF 2, Fig. 6, p. 648
p2p.set_distance(150.0) # km
print 'D, Distance source to site O-C (km): {}'.format(p2p.D)

# z, Zenith angle site, REF 2, Fig. 6, p. 648
p2p.set_zenith(0.0) # deg
print 'z, Site zenith (deg): {}'.format(p2p.get_z_deg())

# beta, Azimuth angle from line of sight to scatter from site, REF 2, Fig. 6, p. 648
p2p.set_beta(0.0) #deg
print 'beta, Relative azimuth line-of-sight to scatter (deg): {}'.format(p2p.get_beta_deg())

# Earth radius of curvature, REF Wikipedia (https://en.wikipedia.org/wiki/Earth_radius)
p2p.set_radius_curvature() #km
print 'R_T, Radius of curvature of the Earth (km): {}'.format(p2p.R_T)

# Scattering distance increment
p2p.set_del_u(0.2) #km
print 'delta_u, Scattering distance increment for finite integration over u (km): {}'.format(p2p.del_u)

## START OF SUM LUMINOUS FLUX FUNCTION
# Unpack point-to-point variables
lat = p2p.lat0
D_OC = p2p.D
zen = p2p.z
beta = p2p.beta
R_T = p2p.R_T
del_u = p2p.del_u

# Set "u", the line of sight to scattering that is integration variable
# u, Direct line distance from scattering point (Q) to observation site (O), REF 2, Fig. 6, p. 648
u_OQ = 25.0 # km
print 'u, Distance scatter to site O-Q (km): {}'.format(u_OQ)

# Constants
# N_m,0 - Molecular density at sea level, REF 2, p. 645
N_m0 = 2.55e19 # cm^-3

# c - Inverse scale altitude, REF. 2, p. 645
c_isa = 0.104 # km^-1

# K - Parameter for relative importance of aerosols to molecules, REF 1, p. 10
K_am = 1.0 # unitless ratio

# sigma_m - Integrated Rayleigh scattering cross-section in the V-band, REF 2, p. 646
sig_m = 1.136e-26 # cm^-2*sr^-1

# Falchi best fit parameters for normalzied emssion
# Weights for angular distributions. REF 1, p. 21
W_a = 1.9e-3
W_b = 5.2e-4
W_c = 7.6e-5
# light loss per night hour after midnight (negative hours before midnight)
d_ll = -0.045

# Equations
print '**EQUATIONS**'
print '*GEOMETRY*'

# Earth angle from source to site, REF 3, p. 308
Chi = D_OC/R_T # radians for trig
Chi_deg = Chi*180/pi # degrees for display
print 'Chi, Arc angle of Earth subtended from source to site OSC (deg): {}'.format(Chi_deg)

# u0, shortest scattering distance based on curvature of the Earth, REF 2, Eq. 21, p. 647
u0 = 2*R_T*sin(Chi/2)**2/(sin(zen)*cos(beta)*sin(Chi)+cos(zen)*cos(Chi)) #km
print 'u0, Shortest scattering distance based on curvature of the Earth (km): {}'.format(u0)

# l, Direct line of sight distance between source and observations site, REF 2, Appendix A (A1), p. 656
l_OC = sqrt(4*R_T**2*sin(Chi/2)**2) # km
print 'l, Direct line of sight distance source to site O-C (km): {}'.format(l_OC)

# q1, Intermediate quantity, REF 2, Appendix A (A1), p. 656, **WITH CORRECTION FROM REF 3, eq. 6, p. 308**
q1 = R_T*(sin(Chi)*sin(zen)*cos(beta) + cos(Chi)*cos(zen) - cos(zen)) # km
print 'q1: {}'.format(q1)

# theta, elevation angle of scatter above source from site (QOC), REF 2, Appendix A (A1), p. 656
theta = arccos(q1/l_OC) # radians
theta_deg = theta*180/pi # degrees
print 'theta, elevation angle of scatter above source from site (QOC) (deg): {}'.format(theta_deg)

## START OF "u" LOOP
# s, Distance from source to scattering CQ, REF 2, Appendix A (A1), p. 656
s_CQ = sqrt((u_OQ - l_OC)**2 + 4*u_OQ*l_OC*sin(Chi/2)**2) # km
print 's, Direct line of sight distance source to scatter O-Q (km): {}'.format(s_CQ)

# h, Height of scattering above Earth reference surface Q, REF 2, Appendix A (A1), p. 656
h_Q = R_T*(sqrt(1 + (u_OQ**2 + 2*u_OQ*R_T*cos(zen))/R_T**2) - 1) # km
print 'h, height of scatter at Q (km): {}'.format(h_Q)

# phi, Elevation angle of emission from C to Q over line of sight from C to O (QCO), REF 2, Appendix A (A1), p. 656
phi = arccos((l_OC**2 + s_CQ**2 - u_OQ**2)/(2*l_OC*s_CQ)) # radians
phi_deg = phi*180/pi # degrees
print 'phi, elevation angle of scatter above site from source (QCO) (deg): {}'.format(phi_deg)

# q3, Intermediate quantity, REF 2, Appendix A (A1), p. 656 
q3 = u_OQ*cos(zen)*cos(Chi) - 2*R_T*sin(Chi/2)**2 # km
print 'q3: {}'.format(q3)

# q2, Intermediate quantity (q2=q3 if zenith angle, z=0), REF 2, Appendix A (A1), p. 656
q2 = u_OQ*sin(zen)*cos(beta)*sin(Chi) + q3 # km
print 'q2: {}'.format(q2)

# Psi, emission angle from source, REF 2, Appendix A (A1), p. 656
Psi = arccos(q2/s_CQ) # radians
Psi_deg = Psi*180/pi # degrees
print 'Psi, emission angle at source (deg): {}'.format(Psi_deg)

# omega, scattering angle at Q, REF 2, Appendix A (A1), p. 656
omega = theta + phi # radians
omega_deg = omega*180/pi # degrees
print 'omega, scattering angle at Q (deg): {}'.format(omega_deg)

print '*ATMOSPHERE*'
# a, Scale height of aerosols, REF 2, p. 646
a_sha = 0.657 + 0.059*K_am # km^-1
print 'a: {}'.format(a_sha)

# p4, Intermediate quantity u-path, REF 2, Appendix A (A2), p. 656
p4 = (a_sha**2*u_OQ**2*cos(zen)**2 + 2*a_sha*u_OQ*cos(zen) + 2)*exp(-a_sha*u_OQ*cos(zen)) - 2
print 'p4: {}'.format(p4)

# p3, Intermediate quantity u-path, REF 2, Appendix A (A2), p. 656
p3 = (c_isa**2*u_OQ**2*cos(zen)**2 + 2*c_isa*u_OQ*cos(zen) + 2)*exp(-c_isa*u_OQ*cos(zen)) - 2
print 'p3: {}'.format(p3)

# p2, Intermediate quantity u-path, REF 2, Appendix A (A2), p. 656
p2 = a_sha**-1*cos(zen*(1 - exp(-a_sha*u_OQ*cos(zen)) + ((16*p4*tan(zen)**2)/(9*pi*2*a_sha*R_T))))**-1
print 'p2: {}'.format(p2)

# p1, Intermediate quantity u-path, REF 2, Appendix A (A2), p. 656
p1 = c_isa**-1*cos(zen*(1 - exp(-c_isa*u_OQ*cos(zen)) + ((16*p3*tan(zen)**2)/(9*pi*2*c_isa*R_T))))**-1
print 'p1: {}'.format(p1)

# ksi1, Extinction of light along u-path from scatter at Q to observation at O, REF 2, Appendix A (A2), p. 656
ksi1 = exp(-N_m0*sig_m*(p1 + 11.778*K_am*p2)) 
print 'ksi1, extinction of light along u-path from scatter at Q to observation at O (unitless): {}'.format(ksi1)

# f4, Intermediate quantity s-path, REF 2, Appendix A (A2), p. 657
f4 = (a_sha**2*s_CQ**2*cos(Psi)**2 + 2*a_sha*s_CQ*cos(Psi) + 2)*exp(-a_sha*s_CQ*cos(Psi)) - 2
print 'f4: {}'.format(f4)

# f3, Intermediate quantity s-path, REF 2, Appendix A (A2), p. 657
f3 = (c_isa**2*s_CQ**2*cos(Psi)**2 + 2*c_isa*s_CQ*cos(Psi) + 2)*exp(-c_isa*s_CQ*cos(Psi)) - 2
print 'f3: {}'.format(f3)

# f2, Intermediate quantity s-path, REF 2, Appendix A (A2), p. 657
f2 = a_sha**-1*cos(Psi*(1 - exp(-a_sha*s_CQ*cos(Psi)) + ((16*f4*tan(Psi)**2)/(9*pi*2*a_sha*R_T))))**-1
print 'f2: {}'.format(f2)

# f1, Intermediate quantity s-path, REF 2, Appendix A (A2), p. 657
f1 = c_isa**-1*cos(Psi*(1 - exp(-c_isa*s_CQ*cos(Psi)) + ((16*f3*tan(Psi)**2)/(9*pi*2*c_isa*R_T))))**-1
print 'f1: {}'.format(f1)

# ksi2, Extinction of light along s-path from emission at C to scatter at Q, REF 2, Appendix A (A2), p. 656
ksi2 = exp(-N_m0*sig_m*(f1 + 11.778*K_am*f2)) 
print 'ksi2, extinction of light along s-path from emission at C to scatter at Q (unitless): {}'.format(ksi2)

# I(Psi), Normalized emission function, MODIFIED FROM REF 1, p. 13 (leaving out natural sky brightness)
I_ne = 1/(2*pi)*(W_a*2*cos(Psi) + W_b*0.554*Psi**4 + W_c*sin(3*Psi))
print 'I(Psi), Normalized emission function at the source: {}'.format(I_ne)

# i(Psi,s), Illuminance per unit flux, REF 2, Eq. 6, p. 644
i_ps = I_ne*ksi2/s_CQ**2
print 'i(Psi,s), Illuminance per unit flux at scattering location Q: {}'.format(i_ps)

# N_m(h), Number density of gaseous component of atmosphere as function of altitude, h, REF 2, Eq. 10, p. 645
N_m = N_m0*exp(-c_isa*h_Q)
print 'N_m(h), Number density of gaseous component of atmosphere as function of altitude, h: {}'.format(N_m)

# N_a*sigma_a, Total integrated scattering cross-section, REF 2, Eq. 12, p. 645 **REARRANGED**
Na_x_siga = K_am*N_m*sig_m*11.11
print 'N_a*sigma_a, Total integrated scattering cross-section of aerosol component as function of altitude, h: {}'.format(Na_x_siga)

# f_m, Angular scattering function for molecular Rayleigh scattering, REF 2, Eq. 13, p. 646
f_m = 3*(1 + cos(omega)**2)/(16*pi)
print 'f_m(omega), Angular scattering function for molecular Rayleigh scattering: {}'.format(f_m)

# f_a, Angular scattering function for aerosol Mie scattering, REF 2, Eq. 14, p. 646
if omega_deg >= 0.0 and omega_deg <= 10.0:
	f_a = 7.5*exp(-0.1249*omega_deg**2/(1 + 0.04996*omega_deg**2))
elif omega_deg > 10.0 and omega_deg <= 124.0:
	f_a = 1.88*exp(-0.07226*omega_deg + 0.0002406*omega_deg**2)
elif omega_deg > 124.0 and omega_deg <= 180.0:
	f_a = 0.025 + 0.015*sin((2.25*omega_deg - 369.0)*(pi/180))
print 'f_a(omega), Angular scattering function for aerosol Mie scattering: {}'.format(f_a)

# S_d, Luminous flux per unit solid angle per unit upward flux (directly from source), REF 2, Eq. 5, p. 644
S_d = (N_m*sig_m*f_m + Na_x_siga*f_a)*i_ps
print 'S_d, Luminous flux per unit solid angle per unit upward flux (directly from source): {}'.format(S_d)

# D_S, Double scattering correction factor, REF 2, Eq. 20, p. 647
D_S = 1 + N_m0*sig_m*(11.11*K_am*f2 + (f1/3))
print 'D_S, Double scattering correction factor: {}'.format(D_S)

# S_u, Total illumance as a function of u, REF 2, Eq. 8, p. 645
S_u = S_d*D_S
print 'S, Total illumance as a function of u: {}'.format(S_u)

# integrand of propogation function, REF 2, Eq. 3, p. 644
df_prop = S_u*ksi1*del_u

# running f_sum = f_sum_previous + df_prop
# exit condition if df_prop < 0.1 percent of f_sum
# print u at the end to check