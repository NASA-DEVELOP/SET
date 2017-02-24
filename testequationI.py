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

# I(Psi)_A = (1/(2*pi))*2*cos(Psi)      # This is the Lambertian circular distribution for stuff reflected from the ground
# I(Psi)_B = (1/(2*pi))*0.554*Psi^4  # This is the horizon lighting model poorly shielded street lighting
# I(Psi)_C = (1/(2*pi))*sin(3*Psi)       # This is the rose petal function that is max at 30 degrees. It is for the intermediate light

#Before solving these equations we need to solve the equation for Psi
# Psi = D/Rt(arcsin*(Rt*sin(D/Rt)/s))

# To solve for Psi we need to solve for s
# s = sqrroot(Rt^2+(Rt+Rs)^2-2*Rt*(Rt+Rs)*cos(D/Rt))
# Where Rt is the average earth radius in kilometers (we can later get more specific to the local average earth radius near our study area)
# Rs = h + Rt where h is the orbital height of NPP from the Earth's surface
# D is the distance from the source to the site

import numpy
from numpy import *

# Input Variables
print '**INPUTS**'

# Latitude of observer
lat_dd = 43.5 # degrees
lat = lat_dd*pi/180 # converted to radians for trig
print 'Latitude (deg): {}'.format(lat_dd)

# Distance from source (C) to observation site (O) along ellipsoid surface
D_OC = 150.0 # km
print 'D, Distance source to site O-C (km): {}'.format(D_OC)

# u, Direct line distance from scattering point (Q) to observation site (O)
u_OQ = 25.0 # km
print 'u, Distance scatter to site O-Q (km): {}'.format(u_OQ)

# z, Zenith angle site
zen_deg = 0.0 # deg
zen = zen_deg*pi/180 # converted to radians for trig
print 'z, Site zenith (deg): {}'.format(zen_deg)

# beta, Azimuth angle from line of sight to scatter from site
beta_deg = 0.0 #deg
beta = beta_deg*pi/180 # converted to radians for trig
print 'beta, Relative azimuth line-of-sight to scatter (deg): {}'.format(beta_deg)

# Constants
# Earth ellipse semi-major orbit
R_equator = 6378.1370 # km (a)

# Earth ellipse semi-minor axis
R_polar = 6356.7523142 # km (b)

# N_m,0 - Molecular density at sea level
N_m0 = 2.55e19 # cm^-3

# c - Inverse scale altitude
c_isa = 0.104 # km^-1

# K - Parameter for relative importance of aerosols to molecules
K_am = 1.0 # unitless ratio

# sigma_m - Integrated Rayleigh scattering cross-section in the V-band
sig_m = 1.136e-26 # cm^-2*sr^-1

# Falchi best fit parameters for normalzied emssion
# Weights for angular distributions
W_a = 1.9e-3
W_b = 5.2e-4
W_c = 7.6e-5
# light loss per night hour after midnight (negative hours before midnight)
d_ll = -0.045

# Equations
print '**EQUATIONS**'
print '*GEOMETRY*'
# Earth radius of curvature
R_T = (R_polar*R_equator**2)/((R_equator*cos(lat))**2+(R_polar*sin(lat))**2)
print 'R_T, Radius of curvature of the Earth (km): {}'.format(R_T)

# Earth angle from source to site
Chi = D_OC/R_T # radians for trig
Chi_deg = Chi*180/pi # degrees for display
print 'Chi, Arc angle of Earth subtended from source to site OSC (deg): {}'.format(Chi_deg)

# u0, shortest scattering distance based on curvature of the Earth
u0 = 2*R_T*sin(Chi/2)**2/(sin(zen)*cos(beta)*sin(Chi)+cos(zen)*cos(Chi)) #km
print 'u0, Shortest scattering distance based on curvature of the Earth (km): {}'.format(u0)

# l, Direct line of sight distance between source and observations site
l_OC = sqrt(4*R_T**2*sin(Chi/2)**2) # km
print 'l, Direct line of sight distance source to site O-C (km): {}'.format(l_OC)

# s, Distance from source to scattering CQ
s_CQ = sqrt((u_OQ - l_OC)**2 + 4*u_OQ*l_OC*sin(Chi/2)**2) # km
print 's, Direct line of sight distance source to scatter O-Q (km): {}'.format(s_CQ)

# h, Height of scattering above Earth reference surface Q
h_Q = R_T*(sqrt(1 + (u_OQ**2 + 2*u_OQ*R_T*cos(zen))/R_T**2) - 1) # km
print 'h, height of scatter at Q (km): {}'.format(h_Q)

# phi, Elevation angle of emission from C to Q over line of sight from C to O (QCO)
phi = arccos((l_OC**2 + s_CQ**2 - u_OQ**2)/(2*l_OC*s_CQ)) # radians
phi_deg = phi*180/pi # degrees
print 'phi, elevation angle of scatter above site from source (QCO) (deg): {}'.format(phi_deg)

# q3, Intermediate quantity 
q3 = u_OQ*cos(zen)*cos(Chi) - 2*R_T*sin(Chi/2)**2 # km
print 'q3: {}'.format(q3)

# q2, Intermediate quantity (q2=q3 if zenith angle, z=0)
q2 = u_OQ*sin(zen)*cos(beta)*sin(Chi) + q3 # km
print 'q2: {}'.format(q2)

# q1, Intermediate quantity
q1 = R_T*(sin(Chi)*sin(zen)*cos(beta) + cos(Chi)*cos(zen) - cos(zen)) # km
print 'q1: {}'.format(q1)

# theta, elevation angle of scatter above source from site (QOC)
theta = arccos(q1/l_OC) # radians
theta_deg = theta*180/pi # degrees
print 'theta, elevation angle of scatter above source from site (QOC) (deg): {}'.format(theta_deg)

# Psi, emission angle from source
Psi = arccos(q2/s_CQ) # radians
Psi_deg = Psi*180/pi # degrees
print 'Psi, emission angle at source (deg): {}'.format(Psi_deg)

# omega, scattering angle at Q
omega = theta + phi # radians
omega_deg = omega*180/pi # degrees
print 'omega, scattering angle at Q (deg): {}'.format(omega_deg)

print '*ATMOSPHERE*'
# a, Scale height of aerosols
a_sha = 0.657 + 0.059*K_am # km^-1
print 'a: {}'.format(a_sha)

# p4, Intermediate quantity u-path
p4 = (a_sha**2*u_OQ**2*cos(zen)**2 + 2*a_sha*u_OQ*cos(zen) + 2)*exp(-a_sha*u_OQ*cos(zen)) - 2
print 'p4: {}'.format(p4)

# p3, Intermediate quantity u-path
p3 = (c_isa**2*u_OQ**2*cos(zen)**2 + 2*c_isa*u_OQ*cos(zen) + 2)*exp(-c_isa*u_OQ*cos(zen)) - 2
print 'p3: {}'.format(p3)

# p2, Intermediate quantity u-path
p2 = a_sha**-1*cos(zen*(1 - exp(-a_sha*u_OQ*cos(zen)) + ((16*p4*tan(zen)**2)/(9*pi*2*a_sha*R_T))))**-1
print 'p2: {}'.format(p2)

# p1, Intermediate quantity u-path
p1 = c_isa**-1*cos(zen*(1 - exp(-c_isa*u_OQ*cos(zen)) + ((16*p3*tan(zen)**2)/(9*pi*2*c_isa*R_T))))**-1
print 'p1: {}'.format(p1)

# ksi1, Extinction of light along u-path from scatter at Q to observation at O
ksi1 = exp(-N_m0*sig_m*(p1 + 11.778*K_am*p2)) 
print 'ksi1, extinction of light along u-path from scatter at Q to observation at O (unitless): {}'.format(ksi1)

# f4, Intermediate quantity s-path
f4 = (a_sha**2*s_CQ**2*cos(Psi)**2 + 2*a_sha*s_CQ*cos(Psi) + 2)*exp(-a_sha*s_CQ*cos(Psi)) - 2
print 'f4: {}'.format(f4)

# f3, Intermediate quantity s-path
f3 = (c_isa**2*s_CQ**2*cos(Psi)**2 + 2*c_isa*s_CQ*cos(Psi) + 2)*exp(-c_isa*s_CQ*cos(Psi)) - 2
print 'f3: {}'.format(f3)

# f2, Intermediate quantity s-path
f2 = a_sha**-1*cos(Psi*(1 - exp(-a_sha*s_CQ*cos(Psi)) + ((16*f4*tan(Psi)**2)/(9*pi*2*a_sha*R_T))))**-1
print 'f2: {}'.format(f2)

# f1, Intermediate quantity s-path
f1 = c_isa**-1*cos(Psi*(1 - exp(-c_isa*s_CQ*cos(Psi)) + ((16*f3*tan(Psi)**2)/(9*pi*2*c_isa*R_T))))**-1
print 'f1: {}'.format(f1)

# ksi2, Extinction of light along s-path from emission at C to scatter at Q
ksi2 = exp(-N_m0*sig_m*(f1 + 11.778*K_am*f2)) 
print 'ksi2, extinction of light along s-path from emission at C to scatter at Q (unitless): {}'.format(ksi2)

# I(Psi), Normalized emission function
I_ne = 1/(2*pi)*(W_a*2*cos(Psi) + W_b*0.554*Psi**4 + W_c*sin(3*Psi))
print 'I(Psi), Normalized emission function at the source: {}'.format(I_ne)

# i(Psi,s), Illuminance per unit flux
i_ps = I_ne*ksi2/s_CQ**2
print 'i(Psi,s), Illuminance per unit flux at scattering location Q: {}'.format(i_ps)

# N_m(h), Number density of gaseous component of atmosphere as function of altitude, h
N_m = N_m0*exp(-c_isa*h_Q)
print 'N_m(h), Number density of gaseous component of atmosphere as function of altitude, h: {}'.format(N_m)

# N_a*sigma_a, Total integrated scattering cross-section
Na_x_siga = K_am*N_m*sig_m*11.11
print 'N_a*sigma_a, Total integrated scattering cross-section of aerosol component as function of altitude, h: {}'.format(Na_x_siga)

# f_m, Angular scattering function for molecular Rayleigh scattering
f_m = 3*(1 + cos(omega)**2)/(16*pi)
print 'f_m(omega), Angular scattering function for molecular Rayleigh scattering: {}'.format(f_m)

# f_a, Angular scattering function for aerosol Mie scattering

