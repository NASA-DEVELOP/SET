"""
Source Name:	Sky Glow Visualization Toolbox.pyt
Author:			Wyoming Cross-Cutting II, NASA DEVELOP National Program
Description:	Python tool to generate artificial sky glow maps using data from NASA
				and NOAA's Visible Infrared Imaging Radiometer Suite (VIIRS) sensor.
"""
###############INCOMPLETE: Testing random functions from Itest.py###############

import arcpy
from numpy import *

# Calculate Gaussian Earth radius of curvature as a function of latitude
def gauss_earth_curvature_radius(center_lat):
	# Earth ellipse semi-major orbit, REF Wikipedia (https://en.wikipedia.org/wiki/Earth_radius)
	R_equator = 6378.1370 # km (a)

	# Earth ellipse semi-minor axis, REF Wikipedia (https://en.wikipedia.org/wiki/Earth_radius)
	R_polar = 6356.7523142 # km (b)

	# Gaussian Earth radius of curvature, REF Wikipedia (https://en.wikipedia.org/wiki/Earth_radius)
	R_curve = ((R_equator**2)*R_polar)/((R_equator*cos(center_lat))**2 + (R_polar*sin(center_lat))**2)

	return R_curve

# Function to streamline parameter creation.
def parameter(dName,name,datatype,defaultValue=None,parameterType=None,direction=None):
	param = arcpy.Parameter(
		displayName = dName,
		name = name,
		datatype = datatype,
		parameterType = "Required",
		direction = "Input")
	param.value = defaultValue
	return param

class Toolbox(object):
	# Define the toolbox.
    def __init__(self):
        self.label = "Sky Glow Visualization Toolbox"
        self.alias = "skyglow"

        # List of tool classes associated with toolbox.
        self.tools = [CreateArtificialSkyGlowMap]

class CreateArtificialSkyGlowMap(object):
	# Define the tool.
	def __init__(self):
		self.label = "Create Artificial Sky Glow Map"
		self.description = """Produce an aritificial sky glow map using data from NASA 
							  and NOAA's Suomi NPP Visible Infrared Imaging Radiometer Suite
							  (VIIRS) Day/Night Band (DNB) sensor."""
		self.canRunInBackground = True

	# Define parameters.
	def getParameterInfo(self):

		# Input VIIRS reference
		vImage = parameter("Input VIIRS Reference","vImage","DEFile")

		# Latitude (Grand Teton National Park = 43.7904 degrees N)
		lat = parameter("Latitude","lat","GPDouble")

		# Length of u for relaxing integration increment (km) ?
		ubr = parameter("Distance from observer to night sky (km)","ubr","GPDouble")

		# Zenith angle
		zen = parameter("Zenith Angle","zen","GPDouble")

		params = [vImage,lat,ubr,zen]
		return params

	# Source code
	def execute(self,params,messages):
		lat = params[0].valueAsText
		ubr = params[1].valueAsText
		zen = params[2].valueAsText
		return
