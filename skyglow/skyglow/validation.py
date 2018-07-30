from osgeo import gdal
from matplotlib import pyplot as plt
from scipy.stats import spearmanr

import os, glob, ntpath

def validate(lat, lon, groundtruth_file, skyglow_folder):
    groundtruth_raster = gdal.Open(groundtruth_file)
    gt_transform = groundtruth_raster.GetGeoTransform()
    gt_data = groundtruth_raster.ReadAsArray()
    gt_x_origin, gt_y_origin = gt_transform[0], gt_transform[3]
    gt_px_width, gt_px_height = gt_transform[1], gt_transform[5]

    gt_vals, vals = [], []
    skyglow_search = os.path.join(skyglow_folder,'*.tif')
    skyglow_tifs = glob.glob(skyglow_search)
    for tif_name in skyglow_tifs:
        tif_name = ntpath.basename(tif_name)
        args_split = tif_name.split('_')
        zenith = float(args_split[3])
        azimuth = float(args_split[4][:-4])
        altitude = abs(zenith - 90)
        print(zenith, azimuth, altitude)
        gt_row = int((altitude - gt_y_origin)/gt_px_height)
        gt_col = int((azimuth - gt_x_origin)/gt_px_width)
        gt_val = gt_data[gt_row][gt_col]
        while gt_val < -1000:
            gt_col = gt_col + 1 if azimuth < 0 else gt_col - 1
            gt_val = gt_data[gt_row][gt_col]
        gt_vals.append(gt_val)

        raster = gdal.Open(os.path.join(skyglow_folder, tif_name))
        transform = raster.GetGeoTransform()
        data = raster.ReadAsArray()
        x_origin, y_origin = transform[0], transform[3]
        px_width, px_height = transform[1], transform[5]
        row = int((lat - y_origin)/px_height)
        col = int((lon - x_origin)/px_width)
        val = data[row][col]
        vals.append(val)
        print(gt_val, val)
    print(gt_vals, vals)
    print(spearmanr(gt_vals, vals))

    fig, ax = plt.subplots()
    ax.set_axisbelow(True)
    ax.grid(linestyle='-', linewidth='0.5', color='gray')
    plt.scatter(x=gt_vals, y=vals, c='r', s=15)
    plt.xlabel('NPS units')
    plt.ylabel('SET units')
    plt.show()

validate(30.31682, -87.26236, 'GroundTruth/anthlightmags_sphere.tif', 'GI_skyglow')
