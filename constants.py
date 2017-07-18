"""
Created: July 17, 2017
Author: Stanley Yu
Description: Module for all constants, links, and text.
"""
SW = 1680
SH = 1050

INSTR = "Welcome to the Skyglow Estimation Toolbox (SET)!\n" \
        "Produce an artificial skyglow map using local inputs and data" \
        "from NASA and NOAA's Suomi NPP Visible Infrared Imaging Radiometer Suite" \
        "(VIIRS) Day/Night Band (DNB) satellite sensor.\n\n" \
        "1. Import VIIRS Image File. Must be a TIFF (.tif) file.\n" \
        "2. Create or import a kernel to convolve images. Input the region's latitude, " \
        "distance at which incrementation speed increases, zenith angle, and azimuth angle to " \
        "create a new kernel. If kernel is already created, input the TIFF file.\n" \
        "3. Click the \"Generate Artificial Skyglow Map\" button. Select destination directory.\n\n" \
        "REF: Garstang, R.H., 1989. Night-sky brightness at observatories and sites. Pub. Astron. Soc. Pac. 101.\n"

ABOUT = "Skyglow Estimation Toolbox (SET) v0.0.1\n" \
        "NASA DEVELOP National Program\n" \
        "===============================================\n\n" \
        "SET Development Contact:\n" \
        "Ryan Avery   ryan.b.avery@nasa.gov\n" \
        "Stanley Yu   stanley98yu@gmail.com\n\n" \
        "WC Summer 2017 Term\n" \
        "Team Lead: Veronica Warda\n" \
        "Ryan Avery, Steven Chao, Stanley Yu\n\n" \
        "WC Spring 2017 Term\n" \
        "Team Leads: Veronica Warda, Benjamin Marcovitz\n" \
        "Aubrey Hilte, Christine Stevens, Eric White\n\n" \
