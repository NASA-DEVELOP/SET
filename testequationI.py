# Our goal is to get a rough calculation of generalized artificial radiance using Dr. Ross's estimation of the light propogation functions, which are sketched out in Cinzano 2000 and Falchi 2016.

# We will first just work on translating these algebraic equations to python code for a given emission angle. We will be using a package called numpy, which has the capability to work with arrays. 
# Arrays are rasters that can't be loaded up into geoprocessing software. They are quick to process and are a simple data format for data that occurs in a matrix.

#Second we will work on integrating each of these algebraic equations to calculate the total amount of radiance occuring from the possible range of emission angles (psi).

#Third, we will apply the integral equations over all "pixels" in the array using a nested for loop (for loop within a for loop).

#Here are our equations

# I(Psi)_A = (1/(2*pi))*2*cos(Psi)      # This is the Lambertian circular distribution for stuff reflected from the ground
# I(Psi)_B = (1/(2*pi))*0.554*Psi^4  # This is the horizon lighting model poorly shielded street lighting
# I(Psi)_C = (1/(2*pi))*sin(3*Psi)       # This is the rose petal function that is max at 30 degrees. It is for the intermediate light
