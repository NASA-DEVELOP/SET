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

def cdm2_to_mag(m):
    """Convert candelas/m2 to mags/arc-second^2.

    Formula taken from http://unihedron.com/projects/darksky/magconv.php
    """
    return log10(m/108000)/(-0.4)

def cdm2_to_nl(m):
    """Convert candelas/m2 to nl.

    Formula taken from https://converter.eu/luminance/#1_Candela/Square_Meter_in_Lambert
    """
    return 0.000314159265358979*m*(1000000000)

def wcm2sr_to_nl(m):
    """Convert W/cm2/sr (NOAA Elvidge data units) to nl.

    Formula taken from https://converter.eu/luminance/#1_Candela/Square_Meter_in_Lambert
    """
    return m*2145.7077824*(1000000000)


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
    gt_x_origin, gt_y_origin = gt_transform[0], gt_transform[3] #coordinates of upper left corner, borders of the pixel
    gt_px_width, gt_px_height = gt_transform[1], gt_transform[5]

    gt_vals, vals, deltas = [], [], []
    skyglow_search = os.path.join(skyglow_folder, '*.tif')
    skyglow_tifs = glob.glob(skyglow_search) #gets all skyglow maps

    #Setting up csv of the data
    download_csv = str(lat) + "val710.csv" #where you want the file to be downloaded to, change name for each run
    csv = open(download_csv, "w")
    columnTitleRow = "azimuth, zenith, lat, long, gt_val, set_val, delta\n"
    csv.write(columnTitleRow)

    #import first into ArcMap and Project Raster (to something like Authalic Sphere).
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
        #This is confusing. Seems like it's checking whether the pixel we're no has an error, and if so going to the next azimuth/pixel.
        #but then wouldn't the azimuth of the checked image no longer correspond to that of the skyglow image we're checking?
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

        #gt_vals.append(gt_val_mcd)
        #this line above was appending nothing, since gt_val_mcd was not a defined function

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
        val_unitconverted = cdm2_to_mag(val_cdm2)  # USE THIS FOR CONVERSION
        #we want to convert to nl to compare with the nl version of NPs data since that is not in log form like mag
        vals.append(val_unitconverted)
        print("val_cdm2, ", val_cdm2, val_unitconverted)

        delta = gt_val - val_unitconverted
        deltas.append(delta)
        print("val_cdm2, val mag, delta:", val_cdm2, val_unitconverted, delta)

        #Adding data to the csv for each skyglow map
        row = str(azimuth) + "," + str(zenith) + "," + str(lat) + "," + str(lon) + "," + str(gt_val) + "," + str(val_unitconverted) + "," + str(delta) + "\n"
        csv.write(row)

    # goodness of fit scores
    #print(spearmanr(gt_vals, vals))
    #print('R-Squared result('+str(r2_score(gt_vals, vals))+')')
    #print('Average delta('+str(np.mean(deltas))+')')

    # draw correlation scatterplot
    fig, ax = plt.subplots()
    ax.set_axisbelow(True)
    ax.grid(linestyle='-', linewidth='0.5', color='gray')
    plt.scatter(x=gt_vals, y=vals, c='r', s=15)
    plt.xlabel('NPS')
    plt.ylabel('SET')
    plt.show()

# Validation on Gulf Islands National Seashore example
validate(38.2198, -78.74026, 'GroundTruth/SHEN_no_trans_mag_sb_gt.tif', 'Park_Skyglow_Maps/SHEN05101418/')
