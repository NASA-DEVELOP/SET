"""
Created: July 17, 2017
Author: Stanley Yu
Description: Module for all constants, links, and text.
"""
SW = 1680
SH = 1050

ICO = './ico/satellite.ico'

INSTR = "Welcome to the Skyglow Estimation Toolbox (SET)!\n" \
        "======================================\n" \
        "Produce an artificial skyglow map using local inputs and data " \
        "from NASA and NOAA's Suomi NPP Visible Infrared Imaging Radiometer Suite" \
        "(VIIRS) Day/Night Band (DNB) satellite sensor.\n\n" \
        "1. Import VIIRS Image File. Must be a TIFF (.tif) file.\n" \
        "2. Create or import a kernel to convolve images. Input the region's latitude, " \
        "distance at which incrementation speed increases, zenith angle, and azimuth angle to " \
        "create a new kernel. If kernel is already created, input the TIFF file.\n" \
        "3. Click the \"Generate Artificial Skyglow Map\" button. Select destination " \
        "directory.\n\n" \
        "Zenith Angle: z (deg)\n" \
        "Azimuth Angle: beta (deg)\n" \
        "Distance at which Integration Speed Increases: u0 (km)\n\n"
        
CDIAG = "REF: Cinzano, P., F. Falchi, C.D. Elvidge and  K.E. Baugh, 2000. The artificial " \
        "night sky brightness mapped from DMSP satellite Operational Linescan System " \
        "measurements. Mon. Not. R. Astron. Soc. 318. Fig. 6\n"

ABOUT = "Skyglow Estimation Toolbox (SET) v0.0.1\n" \
        "NASA DEVELOP National Program\n" \
        "==============================================\n" \
        "SET Development Contact:\n" \
        "Ryan Avery: ryan.b.avery@nasa.gov\n" \
        "Stanley Yu: stanley98yu@gmail.com\n\n" \
        "WC Summer 2017 Term\n" \
        "Team Lead: Veronica Warda\n" \
        "Ryan Avery, Steven Chao, Stanley Yu\n\n" \
        "WC Spring 2017 Term\n" \
        "Team Leads: Veronica Warda, Benjamin Marcovitz\n" \
        "Aubrey Hilte, Christine Stevens, Eric White\n" \
        "==============================================\n\n" \
        "Icon made by Madebyoliver <https://www.flaticon.com/authors/madebyoliver> " \
        "from www.flaticon.com is licensed by CC 3.0 BY <https://creativecommons.org" \
        "/licenses/by/3.0/>" \
