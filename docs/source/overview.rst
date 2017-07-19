=====================
Program Overview
=====================

The Wyoming Cross-Cutting team at NASA DEVELOP aims to create a tool with the ability to accurately model sky glow from artificial sources using publicly available Suomi NPP VIIRS Day/Night Band imagery. We hope that this open source tool informs future research and decision making to mitigate light pollution and to inform the public of the phenomenon.

.. rubric:: Key reasons for the Sky Glow Estimation Toolbox:
* The NPS currently monitors the night sky using in situ surveys with Unihedron Sky Quality Meters. While these cameras take detailed images of the hemisphere, it is costly to travel to these sites and set up the equipment. These observations also only contain information from a single vantage point and at a single point in time. The park service desires a tool that allows them to continually monitor the quality of the night sky at a regional scale.
* Light pollution disrupts circadian rhythms and the production of the hormone melatonin in humans and has been linked to various health disorders such as obesity, tumor growth, depression, and insomnia.
* Light regulates the day-night cycles of plants and animals, therefore light pollution affects the growth of plants, disorients and confuses animals, and impacts breeding cycles.
* With nearly 200,000,000 outdoor light fixtures in the United States, an estimated $2 billion in energy costs is wasted every year due to poorly designed fixtures.


----------------------------------------
Features
----------------------------------------

.. class:: left: VIIRS Day/Night Band Image, right: shapes form around selected areas of the plot

   .. image:: _static/adjusted-convolved-map.png
      :scale: 22%

   .. image:: _static/programShapesActive.png
      :scale: 20%


.. note::
   The Sky Glow Estimation Toolbox is currently in it's early *beta*phase, and requires additional work. See the Github page for a list of Issues we are currently working on.

The Toolbox currently contains a number of great features:

* Calculates skyglow at zenith for a given region
* Saves kernel for a given look angle and region in the form of a GeoTiff
* Calculates skyglow for a given zenith and azimuth angle (this takes more than twice as long as zenith when first run to create the kernel)

A geotiff file from VIIRS Day/Night Band data is read as a numpy array in order to apply the kernel that is created from user-supplied parameters. Users can create kernels for different zenith angles, azimuth angles, regional latitudes, and input VIIRS imagery.

Once a kernel is created, the tool then applies the kernel to the VIIRS data using numpy's Fast Fourier Tranform algorithm to compute the convolved image. Kernel's are stored in their own directory and can be called to make the code run in a matter of seconds rather than hours.

Additional features are still to come!

----------------------------------------
Additional Information
----------------------------------------

* **Study Area:**
     Regional (National Park Scale), can be applied anywhere

* **Earth Observations & Parameters:**
     Suomi NPP VIIRS Day/Night Band

