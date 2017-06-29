"""
Source Name:	Skyglow Estimation Toolbox.pyt
Author:			Wyoming Cross-Cutting II, NASA DEVELOP National Program
Description:	Python tool to generate artificial skyglow maps using data from NASA
				and NOAA's Visible Infrared Imaging Radiometer Suite (VIIRS) sensor.
"""

from arcpy import *
from numpy import *

# Function to streamline parameter creation.
def parameter(dName,name,datatype,paramtype=None,direction=None):
	if paramtype == None:
		paramtype = "Required"
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
        self.tools = [CreateKernel,CreateArtificialSkyglowMap]

class CreateKernel(object):
	# Define the tool
	def __init__(self):
		self.label = "Create Kernel"
		self.description = """Produce a kernel for a specified region. Kernels are required to
							  generate an artificial skyglow map that uses localized parameters
							  such as zenith angle and azimuth angle to depict sources of light
							  pollution."""
		self.canRunInBackground = True

	# Define parameters.
	def getParameterInfo(self):

		# Region Latitude (Grand Teton National Park = 43.7904 degrees N)
		lat = parameter("Region Latitude","lat","GPDouble")

		# Length of u for relaxing integration increment (km) ?
		ubr = parameter("Distance from observer to night sky (km)","ubr","GPDouble")

		# Zenith angle
		zen = parameter("Zenith Angle","zen","GPDouble")

		# Azimuth angle
		azi = parameter("Azimuth Angle","azi","GPDouble")

		# Output kernel path
		kerneltiffpath = parameter("Output Kernel Path","kerneltiffpath","DEFolder",direction="Output")

		params = [lat,ubr,zen,azi,kerneltiffpath]
		return params

	def execute(self,params,messages):
		regionlat_arg = params[0].valueAsText
		ubr_arg = params[1].valueAsText
		zen_arg = params[2].valueAsText
		azimuth_arg = params[3].valueAsText
		kerneltiffpath = params[4].valueAsText

		AddMessage("Program is working.")
		return

class CreateArtificialSkyglowMap(object):
	# Define the tool.
	def __init__(self):
		self.label = "Create Artificial Skyglow Map"
		self.description = """Produce an aritificial sky glow map using data from NASA 
							  and NOAA's Suomi NPP Visible Infrared Imaging Radiometer Suite
							  (VIIRS) Day/Night Band (DNB) sensor."""
		self.canRunInBackground = True

	# Define parameters.
	def getParameterInfo(self):

		# Input VIIRS reference
		vImage = parameter("Input VIIRS Reference","vImage","DEFile")

		params = [vImage]
		return params

	# Source code
	def execute(self,params,messages):
		return