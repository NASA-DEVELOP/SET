# README for User Interface branch#
# as of June 20, 2017
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


### How do I get set up? ###

* Install Anaconda for Python 2 here > https://www.continuum.io/downloads If you already have anaconda for python 3, create an environemnt to work with python 2, see instructions> file:///C:/Users/jvaa/Downloads/conda-cheatsheet.pdf
* Conda comes pre instaled with numpy and many other useful packages for scientific python. The major difference between python 2 and 3 is that print statements are different, which can cause errors.
* That's it! You're ready to run the code. I reccomend Sublime text 3 as a text editor and GitBash as a command prompt
