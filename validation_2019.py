import os
import glob
import ntpath
import gdal

from osgeo import gdal #RL: need permission to install
from matplotlib import pyplot as plt
from scipy.stats import spearmanr
import numpy as np
from math import pi, pow, log10
#from sklearn.metrics import r2_score

# for comparing against NPS data in units of magnitude
def cdm2_to_mag(m):
    """Convert candelas/m2 to mags/arc-second^2.

    Formula taken from http://unihedron.com/projects/darksky/magconv.php
    """
    return log10(m/108000)/(-0.4)

# for comparing against NPS data in units of nanolamberts
def cdm2_to_nl(m):
    """Convert candelas/m2 to nl.

    Formula taken from https://converter.eu/luminance/#1_Candela/Square_Meter_in_Lambert
    """
    return 0.000314159265358979*m*(1000000000)

# DEBUGGING: We were trying to figure out the proper units for SET and checked microcandelas and W/cm2/sr.
# the units are microcandelas, so you shouldn't have to use this other function below
def wcm2sr_to_nl(m):
    """Convert W/cm2/sr (NOAA Elvidge data units) to nl.

    Formula taken from https://converter.eu/luminance/#1_Candela/Square_Meter_in_Lambert
    """
    return m*2145.7077824*(1000000000)

# The main validation function that compares NPS ground-truth to SET skyglow maps (NOT the hemisphere,
# but instead what the hemispheres are built from, ie, hemispheres before the interpolation algorithm)
def validate(lat, lon, groundtruth_file, skyglow_folder):
    """Validate skyglow map data with NPS ground truth hemispherical data.

    lat (float): NPS data latitude
    lon (float): NPS data longitude
    groundtruth_file (str): path to NPS hemisphere data
    skyglow_folder (str): path to skyglow map folder
    """
    groundtruth_raster = gdal.Open(groundtruth_file)
    gt_transform = groundtruth_raster.GetGeoTransform() #get the transformation used to go from pixel units to coordinates
    gt_data = groundtruth_raster.ReadAsArray() #reads raster into memory
    #GetGeoTransform produces an array which is being called here to locate elements of the raster
    #http://www2.geog.ucl.ac.uk/~plewis/geogg122_local/geogg122-old/Chapter4_GDAL/GDAL_Python_bindings.html
    gt_x_origin, gt_y_origin = gt_transform[0], gt_transform[3] #coordinates of upper left corner pixel
    gt_px_width, gt_px_height = gt_transform[1], gt_transform[5] #pixel size in image, ie the step size we'll use to walking through image

    gt_vals, vals, deltas = [], [], []
    skyglow_search = os.path.join(skyglow_folder, '*.tif')
    skyglow_tifs = glob.glob(skyglow_search) #gets all skyglow maps

    #Setting up csv of the data
    download_csv = str(lat) + "_validationOutput.csv" #where you want the file to be downloaded to, change name for each run
    csv = open(download_csv, "w")
    columnTitleRow = "azimuth, zenith, lat, long, gt_val, set_val, delta\n"
    csv.write(columnTitleRow)

    #import first into ArcMap and Project Raster into a spheroid projection (e.g., GCS_Sphere_EMEP)
    #Save that raster as a .tif and that's the input for validation.
    # extract values at lat/lon point for each angle (determined by map file name)
    for tif_name in skyglow_tifs:
        tif_name = ntpath.basename(tif_name)
        args_split = tif_name.split('_')
        zenith = float(args_split[3])
        azimuth = float(args_split[4][:-4])
        altitude = abs(zenith - 90)
        print("zen, azi, alt:", zenith, azimuth, altitude)
        print("y_orig, px_height, x_origin, px_width:", gt_y_origin, gt_px_height, gt_x_origin, gt_px_width)
        print("gt_trans:", gt_transform)

        # get NPS ground truth value: lat-long fixed, we need to find azimuth-zenith corresponding to each skyglow tiff
        # goal is to get the ground-truth values at the areas of the NPS hemisphere we are looking for
        #find which pixel on the ground-truth data corresponds with the skyglow tif you are on
        #Does this work correctly--azimuth altitude not in same units as gt_x_origin, gt_y_origin?
        #guess once I can run this, I can check teh print statement above
        gt_row = int((altitude - gt_y_origin) / gt_px_height)
        print("gt_row:", gt_row)
        if (azimuth >= 180.0):
            azimuth = azimuth - 360
        print("azi:", azimuth)
        gt_col = int((azimuth - gt_x_origin) / gt_px_width)
        print("gt_row, gt_col, altitude:", gt_row, gt_col, altitude)
        print("gt_data shape:", gt_data.shape)
        gt_val = gt_data[gt_row][gt_col]
        print("gt_val", gt_val) #brightness of the corresponding pixel

        # Because the NPS camera equipment may be placed at slightly differently angles each time it is used,
        # the size of the ground-truth image coordinates are sometimes off by a couple pixels
        # this while loop takes care of situations when we may be pulling a null NPS value by the border, by walking us
        # back within the borders of the ground-truth image before saving the pixel value for comparison against SET's value.
        while gt_val < -1000:
            if (90 - altitude) < 0.1:
                gt_row = gt_row + 1
            elif (180 - azimuth) < 0.1:
                gt_col = gt_col -1
            elif abs(-180 - azimuth) < 0.1:
                gt_col = gt_col + 1
            gt_val = gt_data[gt_row][gt_col]

        print("gt_row, col, val:", gt_row, gt_col, gt_val)
        gt_vals.append(gt_val)

        # get SET value for each skyglow tiff
        #We only get the brightness values at the lat-long corresponding to that of the in-situ data
        #We do so for every zenith-azimuth pair, so we end up with values at the red dots in the 3-D hemisphere generated at the ground truth lat-long
        #similar to darksky line 412 area
        raster = gdal.Open(os.path.join(skyglow_folder, tif_name))
        transform = raster.GetGeoTransform()
        data = raster.ReadAsArray()
        x_origin, y_origin = transform[0], transform[3]
        px_width, px_height = transform[1], transform[5]
        row = int((lat - y_origin) / px_height)
        col = int((lon - x_origin) / px_width)
        val = data[row][col]
        val_cdm2 = val*pow(10,-6)  # convert from microcandelas to candelas
        val_unitconverted = cdm2_to_nl(val_cdm2)  # USE THIS FOR CONVERSION; change to cdm2_to_mag if we're using NPS image in units of magnitudes
        #we prefer to convert to nl to compare with the nl version of NPs data since that is not in log form like mag
        vals.append(val_unitconverted)
        print("val_cdm2, ", val_cdm2, val_unitconverted)

        delta = gt_val - val_unitconverted
        deltas.append(delta)
        print("val_cdm2, val mag, delta:", val_cdm2, val_unitconverted, delta)

        #Adding data to the csv for each skyglow map
        row = str(azimuth) + "," + str(zenith) + "," + str(lat) + "," + str(lon) + "," + str(gt_val) + "," + str(val_unitconverted) + "," + str(delta) + "\n"
        csv.write(row)

    # draw correlation scatterplot
    fig, ax = plt.subplots()
    ax.set_axisbelow(True)
    ax.grid(linestyle='-', linewidth='0.5', color='gray')
    plt.scatter(x=vals, y=gt_vals, c='r', s=15)
    plt.xlabel('SET')
    plt.ylabel('NPS')
    plt.show()

# Validation on Gulf Islands National Seashore example
# Note the last two arguments--to keep file paths simple, we put the NPS and SET images we are comparing in subfolders of SET alongside this validation script
validate(38.70301, -109.56808, 'GroundTruth/ARCH_no_trans_nl_sb_gt.tif', 'Park_Skyglow_Maps/ARCH_8.5_mask_0.2u0_ubreak0/')
