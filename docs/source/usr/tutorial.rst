============
**Tutorial**
============

The following tutorial will help guide you on how to use the Skyglow Estimation Toolbox to produce your own artificial skyglow maps. SET's interface is designed to be as easy and simple to use as possible. Don't hesitate to contact us about an issue or suggestion you might have!

**Getting Started**
-------------------

.. image:: ../_static/ui_empty.PNG
   :align: left
   :scale: 45%

Starting the program is a simple two-step process. Making sure that you are in the "skyglow" directory in Git Bash, activate the geospatial environment and then run the program::

  $ source activate geospatial
  $ python Skyglow.py

The image on the left displays the interface that should pop-up when you run the Toolbox. Look at the FAQ or contact us if an issue is preventing you from starting the software. 

In the top left corner of the SET interface are the input options that allow you to choose your VIIRS image and create a kernel using local parameters. The top right corner includes a "Help" button that contains the documentation link, instructions, and "about" window.

|
|

**Getting Data**
----------------

Suomi NPP VIIRS data is distributed to various sources and processed into environmental data records. The National Centers for Environmental Information (NCEI) Earth Observation Group produces monthly average radiance composite images from VIIRS. You can download VIIRS Nighttime Light Data for your region from their website `here`__.

.. note::

   Downloading the data may take a while because the files are so large (up to 1.9GB). Keep this in mind before you hit download!

.. __: https://www.ngdc.noaa.gov/eog/viirs/download_dnb_composites.html

**Importing Files**
-------------------

The first step to creating an artificial skyglow map is to import your VIIRS reference image. SET requires TIFF files in order to perform the proper array calculations. 

.. figure:: ../_static/ui_073117.PNG
   :align: right

   *Image of the SET interface displaying a VIIRS median composite of summer months from 2014-16. The program will use the inputs to generate a kernel and create a map of aritficial sky brightness.*

**Creating a Kernel**
---------------------

Creating a kernel is a matter of defining four arguments: region latitude, atmospheric clarity parameter, zenith angle, and azimuth angle. More information on these parameters can be found in the *Instructions* found in the *Help* menu.

.. note::

    If you have already created a kernel file and wish to use the same parameters as before, simply check the *Import Kernel* option and browse for the file.

**Generating an Artificial Skyglow Map**
----------------------------------------

Once all the necessary parameters are chosen, clicking the "Generate Artificial Skyglow Map" will begin the process of modeling sky brightness. A progress log will appear to help you monitor the generation process. Once the program is finished, a "Finished!" message will display along with a bell sound.
