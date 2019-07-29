"""This module contains functions for estimating light pollution using data from
NASA and NOAA's Visible Infrared Imaging Radiometer Suite (VIIRS) satellite
sensor.

References:
(1) Falchi, F., P. Cinzano, D. Duriscoe, C.C.M. Kyba, C.D. Elvidge, K. Baugh, B.A. Portnov, N.A.
          Rybnikov and R. Furgoni, 2016. The new workd atlas of artificial night sky brightness.
          Sci. Adv. 2.
(2) Cinzano, P., F. Falchi, C.D. Elvidge and  K.E. Baugh, 2000. The artificial night sky
          brightness mapped from DMSP satellite Operational Linescan System measurements. Mon.
          Not. R. Astron. Soc. 318.
(3) Garstang, R.H., 1989. Night-sky brightness at observatories and sites. Pub. Astron. Soc.
          Pac. 101.
"""

from __future__ import division
from __future__ import print_function
import os.path, ntpath
import time
import json
import logging
import sys
import glob
from multiprocessing import Pool
from builtins import range

from numpy import *
from scipy import interpolate
from matplotlib import pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cm
from osgeo import gdal

import skyglow.constants as constants

logger = logging.getLogger()


default_kernel_folder = os.path.join(os.getcwd(), "kernel_lib")
default_output_folder = os.path.join(os.getcwd(), "skyglow_out")

# ISSUES
# check angle fixes ~sgmapper and other
# check beta array fixes ~line 298
# still need to get rid of global variables
# still consider making external variable/function constants file

PsiZ_cond = 89.5*(pi/180)
ubr_arg = 10.0


def main():
    """Main entry point if darksky is run from the command line.

    Calls appropriate function for given action.
    """

    if len(sys.argv) < 2:
        raise ValueError("Usage: python {} <action>. Action can be 'sgmap_single', 'kernel_lib', 'sgmap_multiple', 'hemisphere'".format(sys.argv[0]))

    project_root = os.path.join(os.getcwd(), '/'.join(sys.argv[0].split("/")[:-1]))

    action = sys.argv[1]

    # Estimate a single map of artificial sky brightness with or without
    # an existing 2d propagation function (kernel)
    if action == "sgmap_single":
        try:
            la = float(sys.argv[2])
        except:
            raise ValueError("Latitude argument must be a number")
        try:
            ka = float(sys.argv[3])
        except:
            raise ValueError("Atmospheric clarity argument must be a number")
        try:
            ze = float(sys.argv[4])
        except:
            raise ValueError("Zenith argument must be a number")
        try:
            az = float(sys.argv[5])
        except:
            raise ValueError("Azimuth argument must be a number")
        try:
            fi = os.path.abspath(sys.argv[6])
        except:
            raise ValueError("Specify input VIIRS image")
        sgmapper(la, ka, ze, az, fi)

    # Estimate a group ("library") of 2d propagation functions (kernels)
    elif action == "kernel_lib":
        try:
            la = float(sys.argv[2])
        except:
            raise ValueError("Latitude argument must be a number")
        try:
            ka = float(sys.argv[3])
        except:
            raise ValueError("Atmospheric clarity argument must be a number")
        try:
            fi = os.path.abspath(sys.argv[5])
        except:
            raise ValueError("Specify input VIIRS image")
        # zenith angle list from csv file
        try:
            csv_in = sys.argv[4].split("/")
            csv_path = os.path.join(os.getcwd(), *csv_in)
        # default list of zenith angles from defeault.csv in project root
        except:
            logger.info('CSV file not specified, using default.csv')
            csv_path = project_root + 'default.csv'
        try:
            hem = bool(sys.argv[6])
        except:
            logger.info('Not generating kernels for hemispherical visualization')
            hem = False
        try:
            sync = bool(sys.argv[7])
        except:
            sync = False

        # if asynchronous mode is specified, map kernel generation tasks to CPU cores
        # kernels will be generated concurrently
        if not sync:
            static_args = [(la, ka, fi, hem)]
            with open(csv_path, "rb") as f:
                angle_list = loadtxt(f, delimiter=",",skiprows=1)
            # cartesian product of arguments
            args_product = [(x[0], x[1], x[2], y[0], y[1]) for x in static_args for y in angle_list]
            p = Pool()
            try:
                p.map_async(krn_unpack, args_product).get(9999999)
            except KeyboardInterrupt:
                p.close()
        # synchronous run
        else:
            krnlibber(la, ka, csv_path, fi, hem)

    # Produce a set of artificial skyglow maps based an existing library of propagation functions
    elif action == "sgmap_multiple":
        try:
            fi = os.path.abspath(sys.argv[2])
        except:
            raise ValueError("Specify input VIIRS image")
        # Get kernel library folder path
        try:
            krnfolder = os.abspath(sys.argv[3])
        except:
            logger.info('Kernel folder not specified, using {}'.format(default_kernel_folder))
            krnfolder = default_kernel_folder
        # Get skyglow output folder path
        try:
            outfolder = os.path.abspath(sys.argv[4])
        except:
            logger.info('Output folder not specified, using {}'.format(default_output_folder))
            outfolder = default_output_folder
        multisgmapper(fi, krnfolder, outfolder)

    # Genrate hemisphere from existing skyglow map library
    elif action == "hemisphere":
        try:
            la = float(sys.argv[2])
        except:
            raise ValueError("Latitude argument must be a number")
        try:
            lon = float(sys.argv[3])
        except:
            raise ValueError("Longitude argument must be a number")
        try:
            skyglowfolder = os.abspath(sys.argv[3])
        except:
            raise ValueError("Skyglow map folder must be specified")

        generate_hem(la, lon, skyglowfolder)

    else:
        raise ValueError("No action/incorrect action given. Choose ""sgmap_single"", ""kernel_lib"", or ""sgmap_multiple""")

def sgmapper(centerlat_arg, k_am_arg, zen_arg, azi_arg, filein, prop2filein=""):
    """Create a single artificial brightness map.

    centerlat_arg (float): regional latitude in degrees
    k_am_arg (float): atmospheric clarity ratio
    zen_arg (float): viewing angle from zenith in degrees
    azi_arg (float): viewing angle from North (0 degrees)
    filein (str): VIIRS DNB file path
    prop2filein (str): path to existing kernel tiff (default "")
    """
    kerneltiffpath = prop2filein
    # try to to find existing kernel in the current directory or the kernel_lib directory
    if kerneltiffpath == "":
        # search for default kernel name
        kerneltiffpath = 'kernel_' + str(centerlat_arg) + '_' +  str(k_am_arg) + '_' + str(zen_arg) + '_' + str(azi_arg) + '.tif'
    if (os.path.isfile(kerneltiffpath) is False) and (os.path.isfile(os.path.join("kernel_lib", kerneltiffpath)) is False):
        # If there is no 2d propagation function (kernel), create a new kernel

        # convert angles in degrees to radians
        if abs(centerlat_arg) > 90:
            raise ValueError("Latitude must be between -90 and 90 degrees")
        lat_rad = centerlat_arg*(pi/180)
        azi_rad = azi_arg*(pi/180)
        zen_rad = zen_arg*(pi/180)
        propkernel, totaltime = fsum_2d(lat_rad, k_am_arg, zen_rad, azi_rad, filein)
        kerneltiffpath = 'kernel_' + str(centerlat_arg) + '_' +  str(k_am_arg) + '_' + str(zen_arg) + '_' + str(azi_arg) + '.tif'
        array_to_geotiff(propkernel, kerneltiffpath, filein)
        logger.info("time for prop function ubreak 10: %s", totaltime)
    else:
        krnbase = ntpath.basename(kerneltiffpath)
        sgtags = krnbase.split("_")
        zen_arg = sgtags[3]
        azi_arg = sgtags[4][:-4]
        # Put the skyglow basename together
        # If there is a 2d propagation function (kernel) in the current direcotry and not in the kernel_lib subdirectory,
        # then read it into an array
        if (os.path.isfile(kerneltiffpath) is True) and (os.path.isfile(os.path.join("kernel_lib", kerneltiffpath)) is False):
            kerneldata = gdal.Open(kerneltiffpath)
            propkernel = kerneldata.ReadAsArray()
        # If there is a 2d propagation function (kernel) in the kernel_lib sub directory but not the current directory,
        # then read it into an array
        elif (os.path.isfile(os.path.join("kernel_lib", kerneltiffpath)) == True) and (os.path.isfile(kerneltiffpath) == False):
            kerneldata = gdal.Open(os.path.join("kernel_lib", kerneltiffpath))
            propkernel = kerneldata.ReadAsArray()
        # If there is a 2d propagation function (kernel) in the current directory and subdirectory,
        # then read the one in the current directory into an array
        elif (os.path.isfile(kerneltiffpath) == True) and (os.path.isfile(os.path.join("kernel_lib",kerneltiffpath)) == True):
            kerneldata = gdal.Open(kerneltiffpath)
            propkernel = kerneldata.ReadAsArray()

    # read in VIIRS DNB image
    viirsraster = gdal.Open(filein)
    imagearr = viirsraster.ReadAsArray()

    # DEBUG CODE
    logger.debug('imagearr shape : {}'.format(imagearr.shape))
    logger.debug('imagearr shape rows : {}'.format(imagearr.shape[0]))
    logger.debug('imagearr shape columns : {}'.format(imagearr.shape[1]))
    logger.debug('propkernel shape : {}'.format(propkernel.shape))
    logger.debug('propkernel shape rows : {}'.format(propkernel.shape[0]))
    logger.debug('propkernel shape columns : {}'.format(propkernel.shape[1]))

    # Produce sky glow raster
    skyglowarr = convolve_viirs_to_skyglow(imagearr, propkernel)
    skyglowarr, new_transform = subset_to_200km(skyglowarr, propkernel, filein)
    skyglowpath = ('skyglow_' + ntpath.basename(filein)[:-4] + '_' + str(centerlat_arg) + '_' + str(ubr_arg) + '_'
                   + str(zen_arg) + '_' + str(azi_arg) + 'convolved.tif')
    array_to_geotiff(skyglowarr, skyglowpath, filein, new_transform)
    logger.info("===============\n***Finished!***\n===============\nSkyglow Map saved as:\n" + skyglowpath)
    constants.ding()

def krn_unpack(args):
    """Unpack arguments for parallelized kernel generation."""
    print('Generating kernel', args)
    return generate_krn(*args)

def generate_krn(centerlat_arg, k_am_arg, zenith, azimuth, filein, hem):
    """Generate one kernel for kernel_lib action.

    centerlat_arg (float): regional latitude in degrees
    k_am_arg (float): atmospheric clarity ratio
    zenith (float): viewing angle from zenith in degrees
    azimuth (float): viewing angle from North (0 degrees)
    filein (str): VIIRS DNB file path
    hem (bool): mirror kernel for complementary azimuth angle
    """
    # convert from degrees to radians
    if abs(centerlat_arg) > 90:
        raise ValueError("Latitude must be between -90 and 90 degrees")
    lat_rad = centerlat_arg*(pi/180)
    zen_rad = zenith*(pi/180)
    azi_rad = azimuth*(pi/180)

    outname = 'kernel_{}_{}_{}_{}.tif'
    # if hem is False (don't generate hemispherical kernels) just calculate kernel
    if not hem:
        propkernel, totaltime = fsum_2d(lat_rad, k_am_arg, zenith, azimuth, filein)
    # otherwise find complement angle and flip kernel if complementary kernel file exists
    else:
        # azimuth_complement = -azimuth
        # # complementary kernels do NOT exist for angles 0, 180, and -180
        # if abs(azimuth) != 180.0 and azimuth != 0.0:
        #     azimuth_complement = -azimuth
        #     azimuth_complement_outname = outname.format(centerlat_arg, k_am_arg, zenith, azimuth_complement)
        #     if os.path.isfile(azimuth_complement_outname):
        #         logger.info("Found azimuth complement ({}) kernel, flipping instead of generation".format(azimuth_complement))
        #         complement_kernel = gdal.Open(azimuth_complement_outname)
        #         data = complement_kernel.ReadAsArray()
        #         propkernel, totaltime = fliplr(data), 0
        #     else:
        #         propkernel, totaltime = fsum_2d(lat_rad, k_am_arg, zen_rad, azi_rad, filein)
        # # calculate if complementary does not exist
        # else:
        #     propkernel, totaltime = fsum_2d(lat_rad, k_am_arg, zen_rad, azi_rad, filein)
        propkernel, totaltime = fsum_2d(lat_rad, k_am_arg, zen_rad, azi_rad, filein)
    kerneltiffpath = outname.format(centerlat_arg, k_am_arg, zenith, azimuth)
    array_to_geotiff(propkernel, kerneltiffpath, filein)
    logger.info("time for prop function ({}, {}) ubreak 10: {}".format(zenith, azimuth, totaltime))
    print('Finished kernel', kerneltiffpath)

def krnlibber(centerlat_arg, k_am_arg, angles_file, filein, hem):
    """Create a series of kernels from angle combinations in a csv file.

    centerlat_arg (float): regional latitude in degrees
    k_am_arg (float): atmospheric clarity ratio
    angles_file (str): CSV file with angle combinations to create kernels for
    filein (str): VIIRS DNB file path
    hem (bool): mirror kernel for complementary azimuth angle
    """
    if abs(centerlat_arg) > 90:
        raise ValueError("Latitude must be between -90 and 90 degrees")
    lat_rad = centerlat_arg*(pi/180)
    with open(angles_file, "rb") as f:
        angle_list = loadtxt(f, delimiter=",",skiprows=1)
    for angle_set in angle_list:
        zenith = angle_set[0]
        azimuth = angle_set[1]
        zen_rad = zenith*(pi/180)
        azi_rad = azimuth*(pi/180)

        outname = 'kernel_{}_{}_{}_{}.tif'
        # if hem is False (don't generate hemispherical kernels) just calculate kernel
        if not hem:
            propkernel, totaltime = fsum_2d(lat_rad, k_am_arg, zenith, azimuth, filein)
        # otherwise find complement angle and flip kernel if complementary kernel file exists
        else:
            # azimuth_complement = -azimuth
            # # complementary kernels do NOT exist for angles 0, 180, and -180
            # if abs(azimuth) != 180.0 and azimuth != 0.0:
            #     azimuth_complement = -azimuth
            #     azimuth_complement_outname = outname.format(centerlat_arg, k_am_arg, zenith, azimuth_complement)
            #     if os.path.isfile(azimuth_complement_outname):
            #         logger.info("Found azimuth complement ({}) kernel, flipping instead of generation".format(azimuth_complement))
            #         complement_kernel = gdal.Open(azimuth_complement_outname)
            #         data = complement_kernel.ReadAsArray()
            #         propkernel, totaltime = fliplr(data), 0
            #     else:
            #         propkernel, totaltime = fsum_2d(lat_rad, k_am_arg, zen_rad, azi_rad, filein)
            # # calculate if complementary does not exist
            # else:
            #     propkernel, totaltime = fsum_2d(lat_rad, k_am_arg, zen_rad, azi_rad, filein)
            propkernel, totaltime = fsum_2d(lat_rad, k_am_arg, zen_rad, azi_rad, filein)
        kerneltiffpath = outname.format(centerlat_arg, k_am_arg, zenith, azimuth)
        array_to_geotiff(propkernel, kerneltiffpath, filein)
        logger.info("time for prop function ({}, {}) ubreak 10: {}".format(zenith, azimuth, totaltime))
        print('Finished kernel', kerneltiffpath)

def multisgmapper(filein, krnfolder, outfolder):
    """Create series of artificial brightness maps from kernels in folder.

    filein (str): VIIRS DNB file path
    krnfolder (str): path to folder with kernel tifs
    outfolder (str): location to save skyglow maps to
    """
    # read in VIIRS DNB image
    viirsraster = gdal.Open(filein)
    imgarr = viirsraster.ReadAsArray()

    # Get kernel library file list
    krnsrch = os.path.join(krnfolder,'kernel*.tif')
    krnlist = glob.glob(krnsrch)

    ln = 0
    for krn in krnlist:
        # If the loop has already ran once then rescale imgarr back down to prevent multiple scaling factors being applied
        if ln > 0:
            imgarr *= 10**-9
        ln = 1

        kh = gdal.Open(krn)
        krnarr = kh.ReadAsArray()
        # Break down kernel filename to get elements for skyglow filename
        krnbase = os.path.basename(krn)
        sgtags = krnbase.split("_")
        # Change "kernel" to "skyglow"
        sgtags[0] = "skyglow"
        # Put the skyglow basename together
        sep = "_"
        sgbase = sep.join(sgtags)
        sgfile = os.path.join(outfolder, sgbase)
        # Create skyglow array
        sgarr = convolve_viirs_to_skyglow(imgarr, krnarr)
        sgarr, new_transform = subset_to_200km(sgarr, krnarr, filein)
        # Write skyglow to a geotiff
        array_to_geotiff(sgarr, sgfile, filein, new_transform)

def generate_hem(lat, lon, skyglow_folder):
    """Generate skyglow hemisphere for observation point within skyglow map.

    Hemisphere is generated from multiple skyglow maps at same region, each map
    being for a different viewing angle. The value for each angle is taken from
    a skyglow map and the hemisphere is built by interpolating between the
    extracted values. There is a minimum number of values (hence maps) that need
    to exist for proper interpolation.

    lat (float): observation point latitude
    lon (float): observation point longitude
    skyglow_folder (str): path to folder with skyglow maps
    """
    if abs(lat) > 90:
        raise ValueError("Latitude must be between -90 and 90 degrees")
    zen, azi, vals = [], [], []
    # list all tiff files in skyglow_folder (that start with skyglow)
    skyglow_search = os.path.join(skyglow_folder,'skyglow*.tif')
    skyglow_tifs = glob.glob(skyglow_search)
    for tif_name in skyglow_tifs:
        # extract zenith/azimuth from tif name
        tif_name = ntpath.basename(tif_name)
        args_split = tif_name.split('_')
        zenith = float(args_split[3])
        azimuth = float(args_split[4][:-4])
        zen.append(zenith)
        azi.append(azimuth)

        # open tif, read as array, and extract pixel value
        raster = gdal.Open(os.path.join(skyglow_folder, tif_name))
        transform = raster.GetGeoTransform()
        data = raster.ReadAsArray()
        x_origin, y_origin = transform[0], transform[3]
        px_width, px_height = transform[1], transform[5]
        row = int((lat - y_origin)/px_height)
        col = int((lon - x_origin)/px_width)
        try:
            val = data[row][col]
        except:
            raise ValueError("Latitude/longitude must be within bounds of skyglow map.")
        # remove anomalous values near azimuth 180 (why is it happening?)
        #if azimuth == 180 or azimuth == -180:
        #    val = 0
        vals.append(val)
        print(zenith, azimuth, val)

    # I chose Rbf as the interpolator because of https://stackoverflow.com/a/37872172
    interpolator = interpolate.Rbf(azi, zen, vals, function='cubic')
    azinew = arange(-180, 181, 1)
    zennew = arange(0, 81, 1)
    azi_dense, zen_dense = meshgrid(linspace(-180, 180, 360), linspace(0, 80, 80))
    z_dense = interpolator(azi_dense, zen_dense)

    # convert zenith,azimuth, and values in square grid to Hammer projection
    lat = abs(subtract(zen, 90))
    x_hammer, y_hammer = to_hammer_xy(lat, azi)
    y_hammer = subtract(abs(subtract(y_hammer, 90)), 10)
    z_hammer = to_hammer_z(z_dense)

    # fill in nan values within hemisphere (but only inside) if any exist
    for row in range(z_hammer.shape[0]):
        ind = where(~isnan(z_hammer[row]))[0]
        if ind.size != 0:
            first, last = ind[0], ind[-1]
            mask = isnan(z_hammer[row,:])
            mask[:first] = 0
            mask[last:] = 0
            z_hammer[row, mask] = interp(flatnonzero(mask), flatnonzero(~mask), z_hammer[row, ~mask])

    # draw the resulting figure
    # UNCOMMENT BELOW TO USE NPS CUSTOM HEMISPHERE COLORMAP (ALSO DELETE CMAP VARIABLE)
    #RGB values originally from magnitudes.lyr
    # colormap_file = 'colormap_magnitudeslyr.txt'
    # mag_start, mag_end, R, G, B = np.loadtxt(colormap_file).T    #RGB in 0-255 scale
    #
    # #color positions in 0-1 scale
    # mag_range = mag_start[-1]-mag_start[0]
    # mag_percent = (mag_start-mag_start[0])/mag_range
    #
    # #RGB values in 0-1 scale
    # R /= 255.
    # G /= 255.
    # B /= 255.
    #
    # #color map lists
    # red = []; green = []; blue = []
    # for i in range(len(mag_start)):
    #     red.append((mag_percent[i],R[i],R[i]))
    #     green.append((mag_percent[i],G[i],G[i]))
    #     blue.append((mag_percent[i],B[i],B[i]))
    #
    # #declare color map setting
    # cdict = {'red':red,'green':green,'blue':blue}
    #
    # #register the color map
    # plt.register_cmap(name='NPS_mag', data=cdict)

    fig = plt.figure(1, figsize=(10, 6), dpi=100)
    ax = fig.add_subplot(111)
    ax.set_facecolor('white')
    cmap = cm.rainbow
    cmap.set_bad('white',1)
    cax = ax.imshow(z_hammer, extent=(-180, 180, 80, 0), interpolation='nearest', cmap=cmap)
    cbar = fig.colorbar(cax, ticks=[nanmin(z_hammer), (nanmax(z_hammer)-nanmin(z_hammer))/2, nanmax(z_hammer)], orientation='horizontal')
    cbar.set_label('Sky brightness')
    cbar.ax.set_xticklabels(['Low', 'Medium', 'High'])  # horizontal colorbar
    ax.scatter(x=x_hammer, y=y_hammer, c='r', s=10) #the red dots marking out sample points
    plt.show()

def to_hammer_xy(lat, lon):
    """Convert geographic point to Hammer projection.

    Based directly on formulas in https://en.wikipedia.org/wiki/Hammer_projection

    lat (float): latitude of point
    lon (float): longitude of point
    """
    lat_rad = abs(multiply(lat, pi/180))
    lon_rad = multiply(lon, pi/180)
    denominator = sqrt(1+cos(lat_rad)*cos(lon_rad/2))
    x = (2*sqrt(2)*cos(lat_rad)*sin(lon_rad/2))/denominator
    y = (sqrt(2)*sin(lat_rad))/denominator
    return rint(multiply(x, 180/pi)), rint(multiply(y, 180/pi))

def to_hammer_z(values):
    """Convert 2d grid (matrix) to Hammer projection.

    Moves the values within the matrix to make a hemisphere shape.
    Based directly on formulas in https://en.wikipedia.org/wiki/Hammer_projection

    values (array): 2D numpy array of values
    """
    vals_shape = values.shape
    new_values = empty(vals_shape)  # create array to fill with values in new places
    new_values[:] = nan
    for i in range(vals_shape[0]):
        lat_rad = abs(i-90)*pi/180
        for j in range(vals_shape[1]):
            lon_rad = (j-180)*pi/180

            # project from lat/lon to x/y in Hammer
            denominator = sqrt(1+cos(lat_rad)*cos(lon_rad/2))
            x = (2*sqrt(2)*cos(lat_rad)*sin(lon_rad/2))/denominator
            y = (sqrt(2)*sin(lat_rad))/denominator

            # convert back from rads to degrees
            x = int(round(x*180/pi)) + 180
            y = int(abs(round(y*180/pi) - 90))

            # set value in new array in appropriate place
            new_values[y,x] = values[i,j]
    return new_values

def convolve_viirs_to_skyglow(dnbarr, proparr):
    """Convolve VIIRS DNB with light propagation kernel.

    Convolution happens via multiplication in Fourier space (frequency domain).
    An explanation: https://en.wikipedia.org/wiki/Convolution_theorem

    dnbarr (array): 2D numpy array of VIIRS DNB data
    proparr (array): 2D light propagation array
    """
    # Convert to float 32 for Fourier, scale, and round
    # falchi assumed natural sky brightness to be 174 micro cd/m^2 = 2.547e-11 watt/cm^2/steradian (at 555nm)
    # Not sure if this is correct scaling factor, I assume that this makes the output prop image in units of cd/m^2
    viirs_scaling_factor = 10**9
    dnbarr *= viirs_scaling_factor
    proparr = float32(nan_to_num(proparr))

    # DEBUG CODE
    logger.debug('dnbarr.shape[0] : {}'.format(dnbarr.shape[0]))
    logger.debug('dnbarr.shape[1] : {}'.format(dnbarr.shape[1]))
    logger.debug('proparr.shape[0] : {}'.format(proparr.shape[0]))
    logger.debug('proparr.shape[1] : {}'.format(proparr.shape[1]))

    # 400km buffer check using kernel size which is supposed to be 400km radius
    if dnbarr.shape[0] < proparr.shape[0] or dnbarr.shape[1] < proparr.shape[1]:
        raise ValueError("VIIRS image is too small for kernel. Study area must have 400km buffer around it.")

    # generalized padding of kernel so that fft can run
    pad_left = (dnbarr.shape[0] - proparr.shape[0])//2
    pad_right = (dnbarr.shape[0] - proparr.shape[0])//2
    if dnbarr.shape[0] % 2 == 0:
        pad_right += 1
    pad_up = (dnbarr.shape[1] - proparr.shape[1])//2
    pad_down = (dnbarr.shape[1] - proparr.shape[1])//2
    if dnbarr.shape[1] % 2 == 0:
        pad_down += 1
    padded_prop = pad(proparr,((pad_left,pad_right),(pad_up,pad_down)), 'constant', constant_values = 0)

    # DEBUG CODE
    logger.debug('padded_prop : {}'.format(array(padded_prop).shape))
    logger.debug('padded_prop[0] : {}'.format(array(padded_prop).shape[0]))
    logger.debug('padded_prop[1] : {}'.format(array(padded_prop).shape[1]))
    logger.debug('pad_left : {}'.format(array(pad_left)))
    logger.debug('pad_right : {}'.format(array(pad_right)))
    logger.debug('pad_up : {}'.format(array(pad_up)))
    logger.debug('pad_down : {}'.format(array(pad_down)))

    # propagation function applied via fft must be rotated 180 degrees
    padded_prop = fliplr(flipud(padded_prop))

    # for convolution FFT comparison
    # subset = 50
    # prows = padded_prop.shape[0]
    # pcols = padded_prop.shape[1]
    # irows = padded_prop.shape[0]
    # icols = padded_prop.shape[1]

    # padded_prop = padded_prop[(prows//2)-subset:(prows//2)+subset, (pcols//2)-subset:(pcols//2)+subset]
    # imagearr = imagearr[(irows//2)-subset:(irows//2)+subset, (icols//2)-subset:(icols//2)+subset]

    # Fourier Transform Method
    np_dft_prop_im = fft.fft2(padded_prop)
    np_dft_kernel_shift = fft.fftshift(np_dft_prop_im)
    np_magnitude_spectrum = 20*log(abs(np_dft_kernel_shift))

    # DEBUG CODE
    logger.debug('np_dft_prop_im shape : {}'.format(np_dft_prop_im.shape))
    logger.debug('np_dft_kernel_shift shape : {}'.format(np_dft_kernel_shift.shape))
    logger.debug('np_magnitude_spectrum shape : {}'.format(np_magnitude_spectrum.shape))

    np_dft_viirs_im = fft.fft2(dnbarr)
    np_dft_viirs_shift = fft.fftshift(np_dft_viirs_im)
    np_magnitude_spectrum_viirs = 20*log(abs(np_dft_viirs_shift))

    # ADDED DEBUG CODE
    logger.debug('np_dft_viirs_im shape : {}'.format(np_dft_viirs_im.shape))
    logger.debug('np_dft_viirs_shift shape : {}'.format(np_dft_viirs_shift.shape))
    logger.debug('np_magnitude_spectrum_viirs shape : {}'.format(np_magnitude_spectrum_viirs.shape))

    kernel_inv_shift = fft.ifftshift(np_dft_kernel_shift)
    viirs_inv_shift = fft.ifftshift(np_dft_viirs_shift)

    # ADDED DEBUG CODE
    logger.debug('kernel_inv_shift shape : {}'.format(kernel_inv_shift.shape))
    logger.debug('viirs_inv_shift shape : {}'.format(viirs_inv_shift.shape))

    FFT_product_inverse = abs(fft.fftshift(fft.ifft2(kernel_inv_shift * viirs_inv_shift)))
    return FFT_product_inverse

    # Comparison with Slow Convolution (Make sure to subset first) these give slightly different answers
    # apply kernel: https://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.ndimage.filters.convolve.html
    # filtered = ndimage.convolve(imagearr, padded_prop, mode='constant', cval = 0.0)
    # plt.subplot(121),plt.imshow(filtered, norm = colors.LogNorm(), cmap = 'gray')
    # plt.title('Slow Convolution Product'), plt.xticks([]), plt.yticks([])
    # plt.subplot(122),plt.imshow(FFT_product_inverse, norm = colors.LogNorm(), cmap = 'gray')
    # plt.title('Fast FFT Product'), plt.xticks([]), plt.yticks([])
    # plt.show()
    # ##############################################################################

def subset_to_200km(convolved_img, kernel, original_viirs_image):
    """Remove 200km from border of skyglow map.

    convolved_img (array): skyglow map
    kernel (array): light propagation kernel
    original_viirs_image (array): VIIRS DNB data
    """
    # 200km north to south vs east to west may span different number of cells
    km200ud = int(ceil(kernel.shape[1]/2))
    km200lr = int(ceil(kernel.shape[0]/2))
    subsetted = convolved_img[km200ud:-km200ud,km200lr:-km200lr]
    # make sure to move x and y origin to subsetted corner
    imdata = gdal.Open(original_viirs_image)
    transform = list(imdata.GetGeoTransform())
    transform[0] += km200lr*transform[1]
    transform[3] += km200ud*transform[5]
    return subsetted, transform

def fsum_2d(cenlat, k_am, zen, azi, fin):
    """Create 2D light propagation function.

    cenlat (float): latitude in radians
    k_am (float): atmospheric clarity ratio
    zen (float): zenith angle in radians
    azi (float): zaimuth angle in radians
    fin (str): VIIRS DNB file path
    """
    # Input Variables
    logger.info('**INPUTS**')

    # 2d propagation function pixel size in degrees and radians
    p_deg = .0041666667
    p_rad = p_deg*pi/180.0

    # args in degree form
    cenlat_deg = cenlat*(180/pi)
    zen_deg = zen*(180/pi)
    azi_deg = azi*(180/pi)

    # central latitude for estimation of the 2d propagation function
    logger.info('Central latitude (deg): {}'.format(cenlat*(180/pi)))

    # central latitude for estimation of the 2d propagation function
    logger.info('Atmospheric clarity indicator: {}'.format(k_am))

    # z, Zenith angle site, REF 2, Fig. 6, p. 648
    logger.info('z, Site zenith (deg): {}'.format(zen*(180/pi)))

    # Azimuth angle site, REF 2, Fig 6 p. 648
    logger.info('Site azimuth as referenced from North (deg): {}'.format(azi*(180/pi)))

    # ubr, Length of u for relaxing integration increment
    ubr = ubr_arg
    logger.info('ubr, Length of u for relaxing integration increment (km): {}'.format(ubr))

    # Gaussian Earth radius of curvature, REF Wikipedia (https://en.wikipedia.org/wiki/Earth_radius)
    R_T = gauss_earth_curvature_radius(cenlat)
    logger.info('R_T, Radius of curvature of the Earth (km): {}'.format(R_T))

    # Create latitude/longitude arrays
    source_lat, rltv_long, cent_row, cent_col = create_latlon_arrays(R_T, cenlat, p_rad)

    # DEBUGGING: Printing out latitude and longitude arrays for each park. For each park lat and long arrays should be the same across zenith and azimuths.
    #logger.info('file path: {}'.format(fin))
    #source_lat_path = 'source_lat_asymmetrytest_' + str(cenlat_deg) + '_' +  str(zen) + '_' + str(azi) + '.tif'
    # array_to_geotiff(source_lat, source_lat_path, fin)
    # rltv_long_path = 'rltv_long_asymmetrytest_' + str(cenlat_deg) + '_' +  str(zen) + '_' + str(azi) + '.tif'
    # array_to_geotiff(rltv_long, rltv_long_path, fin)

    # Distance from source (C) to observation site (O) along ellipsoid surface, REF 2, Fig. 6, p. 648
    # using haversine formula
    D_OC = 2.0*R_T*arcsin(sqrt(sin((source_lat - cenlat)/2.0)**2.0 + cos(cenlat)*cos(source_lat)*sin(rltv_long/2.0)**2.0))

    # Assignment of NaNs or null values outside of 200 km
    D_OC[D_OC > 201.0] = NaN
    # Delete NaN rows and columns
    D_OC = D_OC[~isnan(D_OC).all(axis=1)]
    D_OC = D_OC[:, ~isnan(D_OC).all(axis=0)]

    # Check of D_OC shape after assigning NaNs outside of 200 km radius
    logger.info('kernel width in pixels, trimmed: {}'.format(D_OC.shape[1]))
    logger.info('kernel height in pixels, trimmed: {}'.format(D_OC.shape[0]))
    widthcenter = (D_OC.shape[1] + 1)//2
    heigthcenter = (D_OC.shape[0] + 1)//2
    # reassignment of center value, need to use better method
    D_OC[cent_row,cent_col] = .01

    # Earth angle from source to site, REF 3, p. 308\
    Chi = D_OC/R_T

    # beta angle from source to site
    beta = create_beta(source_lat, rltv_long, Chi, azi, cent_row, cent_col)

    # u0, shortest scattering distance based on curvature of the Earth, REF 2, Eq. 21, p. 647
    u0 = 2.0*R_T*sin(Chi/2.0)**2.0/(sin(zen)*cos(beta)*sin(Chi)+cos(zen)*cos(Chi)) #km
    # DEBUGGING: Example of how to print a variable to geotiff format for debugging inspection
    #array_to_geotiff(u0, "u0asymmetrytest.tif", fin)

    # l, Direct line of sight distance between source and observations site, REF 2, Appendix A (A1), p. 656
    # l_OC and D_OC are similar as expected
    l_OC = sqrt(4.0*(R_T**2.0)*((sin(Chi/2.0))**2.0)) # km

    # q1, Intermediate quantity, REF 2, Appendix A (A1), p. 656, **WITH CORRECTION FROM REF 3, eq. 6, p. 308**
    q1 = R_T*(sin(Chi)*sin(zen)*cos(beta) + cos(Chi)*cos(zen) - cos(zen)) # km

    # theta, elevation angle of scatter above source from site (QOC), REF 2, Appendix A (A1), p. 656
    theta = arccos(q1/l_OC) # radians

    # simulates what a 180 degrees azimuth input would look like to the program to make sure that rounding errors dont prevent the if statement
    logger.info('azi, {}'.format(azi))
    azi_of_180 = float(180)*pi/180.0
    #DEBUGGING: Commented out mirroring in case it was causing problem at 180 azimuth
    # if azi in {0.0, azi_of_180, -azi_of_180}:
    #     # Get left arrays to cut processing time in half
    #     Chileft = Chi[0:, 0:widthcenter]
    #     u0left = u0[0:, 0:widthcenter]
    #     l_OCleft = l_OC[0:, 0:widthcenter]
    #     thetaleft = theta[0:, 0:widthcenter]
    #     betaleft = beta[0:, 0:widthcenter]
    #     PropSumArrayleft = zeros_like(l_OCleft)
    #     kerdim = l_OCleft.shape
    #     ker10per = kerdim[0]//10
    #
    #     logger.info("Beginning kernel...")
    #     start = time.time()
    #
    #     # 2d iteration for integrating from u0 to infinity to create propagation function for each element
    #     #for p,c,u,l,t in itertools.izip(nditer(PropSumArrayleft, op_flags=['readwrite']),nditer(Chileft, op_flags=['readwrite']),nditer(u0left, op_flags=['readwrite']), nditer(l_OCleft, op_flags=['readwrite']),nditer(thetaleft, op_flags=['readwrite'])):
    #     #    p[...] = fsum_single(R_T, c, u, l, t, 0.0, zen, ubr, k_am)
    #
    #     for ii in range(kerdim[0]):
    #         if mod(ii, ker10per) == 0:
    #             interm = time.time() - start
    #             logger.info("Kernel (%d, %d), %d percent complete, %d minutes",
    #                         zen, azi, ceil(100.0*ii/kerdim[0]), interm/60.0)
    #         for jj in range(kerdim[1]):
    #             PropSumArrayleft[ii][jj], return2 = fsum_single(R_T, Chileft[ii][jj],
    #                                                    u0left[ii][jj],
    #                                                    l_OCleft[ii][jj],
    #                                                    thetaleft[ii][jj],
    #                                                    betaleft[ii][jj],
    #                                                    zen, ubr, k_am)
    #             if (ii == 0 and jj==0):
    #                 print("return1: " + str(PropSumArrayleft[ii][jj]) + "return2: " + str(return2))
    #
    #     time_kern = time.time() - start
    #     logger.info("Time to produce kernel: %d minutes", time_kern/60.0)
    #
    #     # Mirror left side of kernel array to the right
    #     PropSumArrayright = fliplr(PropSumArrayleft[:, 1:])
    #
    #     # Complete 2d propagation function by putting left and right together
    #     PropSumArray = hstack((PropSumArrayleft, PropSumArrayright))
    # else:
    #     # debug checks for subset arrays
    #     logger.debug('Problem Index Chi: %s', Chi[523, 470])
    #     logger.debug('Chi: %s', Chi)
    #     logger.debug('Problem Index u0: %s', u0[523, 470])
    #     logger.debug('u0: %s', u0)
    #     logger.debug('Problem Index l_OC: %s', l_OC[523, 470])
    #     logger.debug('l_OC: %s', l_OC)
    #     logger.debug('Problem Index theta: %s', theta[523, 470])
    #     logger.debug('theta: %s', theta)
    #     logger.debug('Surrounding Indices beta: %s', beta[522:525, 470])
    #     logger.debug('Problem Index beta: %s', beta[523, 470])
    #     logger.debug('beta: %s', beta)

    # initialize kernel matrix for 2d iteration
    PropSumArray = zeros_like(l_OC)

    # DEBUGGING: Initializing matrices for debugging fsum_single outputs which feed into fsum_2d
    Su_max_array = zeros_like(l_OC)
    Su_min_array = zeros_like(l_OC)
    Su_u0val_array = zeros_like(l_OC)

    Sd_max_array = zeros_like(l_OC)
    Sd_min_array = zeros_like(l_OC)
    Sd_u0val_array = zeros_like(l_OC)

    DS_max_array = zeros_like(l_OC)
    DS_min_array = zeros_like(l_OC)
    DS_u0val_array = zeros_like(l_OC)

    ips_max_array = zeros_like(l_OC)
    ips_min_array = zeros_like(l_OC)
    ips_u0val_array = zeros_like(l_OC)

    kerdim = l_OC.shape
    ker10per = kerdim[0]//10

    logger.info("Beginning kernel...")
    start = time.time()

    # 2d iteration for integrating from u0 to infinity to create propagation function for each element
    #for p,c,u,l,t,b in itertools.izip(nditer(PropSumArray, op_flags=['readwrite']), nditer(Chi, op_flags=['readwrite']), nditer(u0, op_flags=['readwrite']), nditer(l_OC, op_flags=['readwrite']), nditer(theta, op_flags=['readwrite']), nditer(beta, op_flags=['readwrite'])):
    #    p[...] = fsum_single(R_T, c, u, l, t, b, zen, ubr, k_am)

    for ii in range(kerdim[0]):
        if mod(ii, ker10per) == 0:
            interm = time.time() - start
            logger.info("Kernel (%d, %d), %d percent complete, %d minutes",
                        zen, azi, ceil(100.0*ii/kerdim[0]), interm/60.0)
        for jj in range(kerdim[1]):
            # returning a dictionary with value for this pixel of the kernel plus debugging values
            results_packed = fsum_single(R_T, Chi[ii][jj], u0[ii][jj],
                                               l_OC[ii][jj], theta[ii][jj],
                                               beta[ii][jj], zen, ubr, k_am)
            # unpacking results dictionary to separate kernel values from debugging values
            # if statements address cases where the kernel value is null due to l_OC being null
            if abs(results_packed["total_sum"] - 88888) < 0.1:
                # assigning kernel value to kernel output matrix
                PropSumArray[ii][jj] = NaN
                # DEBUGGING: remainder of this if statement is assigning values to variables for the debugging of fsum_single
                Su_max_array[ii][jj] = NaN
                Su_min_array[ii][jj] = NaN
                Su_u0val_array[ii][jj] = NaN

                Sd_max_array[ii][jj] = NaN
                Sd_min_array[ii][jj] = NaN
                Sd_u0val_array[ii][jj] = NaN

                DS_max_array[ii][jj] = NaN
                DS_min_array[ii][jj] = NaN
                DS_u0val_array[ii][jj] = NaN

                ips_max_array[ii][jj] = NaN
                ips_min_array[ii][jj] = NaN
                ips_u0val_array[ii][jj] = NaN

            else:
                # assigning kernel value at the pixel to the kernel output matrix for non-NaN values
                PropSumArray[ii][jj] = results_packed["total_sum"]
                # DEBUGGING: remainder of this if statement is assigning values to variables for the debugging of fsum_single
                Su_max_array[ii][jj] = results_packed["Su_max"]
                Su_min_array[ii][jj] = results_packed["Su_min"]
                Su_u0val_array[ii][jj] = results_packed["Su_u0_val"]

                Sd_max_array[ii][jj] = results_packed["Sd_max"]
                Sd_min_array[ii][jj] = results_packed["Sd_min"]
                Sd_u0val_array[ii][jj] = results_packed["Sd_u0_val"]

                DS_max_array[ii][jj] = results_packed["DS_max"]
                DS_min_array[ii][jj] = results_packed["DS_min"]
                DS_u0val_array[ii][jj] = results_packed["DS_u0_val"]

                ips_max_array[ii][jj] = results_packed["ips_max"]
                ips_min_array[ii][jj] = results_packed["ips_min"]
                ips_u0val_array[ii][jj] = results_packed["ips_u0_val"]
            # log print to quickly check kernel value at the center pixel and debugging matrix values
            if (ii == kerdim[0]/2 and jj== kerdim[1]/2):
                # printing out kernel value at center pixel as a quick check
                print("PropSumArray_val: " + str(PropSumArray[ii][jj]))
                #DEBUGGING: printing out debugging values at center for reference
                print("Su_max: " + str(Su_max_array[ii][jj]))
                print("Su_min: " + str(Su_min_array[ii][jj]))
                print("Sd_max: " + str(Sd_max_array[ii][jj]))
                print("Sd_min: " + str(Sd_min_array[ii][jj]))
                print("DS_max: " + str(DS_max_array[ii][jj]))
                print("DS_min: " + + str(DS_min_array[ii][jj]))
                print("ips_max: " + str(ips_max_array[ii][jj]))
                print("ips_min: " + + str(ips_min_array[ii][jj]))

    time_kern = time.time() - start
    logger.info("Time to produce kernel: %d minutes", time_kern/60.0)

    logger.debug('Surrounding Indices PropSumArray: %s', PropSumArray[522:525, 470])
    logger.debug('Problem Index PropSumArray: %s', PropSumArray[523, 470])
    logger.debug('PropSumArray: %s', PropSumArray)

    # DEBUGGING: printing out geotiffs of debugging matrices with maximum, minimum, and integration
    # starting value at each pixel for each potentially problem fsum_single variable
    Su_max_path = 'Su_max' + str(cenlat_deg) + '_' +  str(zen_deg) + '_' + str(azi_deg) + '.tif'
    array_to_geotiff(Su_max_array, Su_max_path, fin)
    Su_min_path = 'Su_min' + str(cenlat_deg) + '_' +  str(zen_deg) + '_' + str(azi_deg) + '.tif'
    array_to_geotiff(Su_min_array, Su_min_path, fin)
    Su_u0val_path = 'Su_u0val' + str(cenlat_deg) + '_' +  str(zen_deg) + '_' + str(azi_deg) + '.tif'
    array_to_geotiff(Su_u0val_array, Su_u0val_path, fin)
    print("Created Su geotiffs")

    Sd_max_path = 'Sd_max' + str(cenlat_deg) + '_' +  str(zen_deg) + '_' + str(azi_deg) + '.tif'
    array_to_geotiff(Sd_max_array, Sd_max_path, fin)
    Sd_min_path = 'Sd_min' + str(cenlat_deg) + '_' +  str(zen_deg) + '_' + str(azi_deg) + '.tif'
    array_to_geotiff(Sd_min_array, Sd_min_path, fin)
    Sd_u0val_path = 'Sd_u0val' + str(cenlat_deg) + '_' +  str(zen_deg) + '_' + str(azi_deg) + '.tif'
    array_to_geotiff(Sd_u0val_array, Sd_u0val_path, fin)
    print("Created Sd geotiffs")

    DS_max_path = 'DS_max' + str(cenlat_deg) + '_' +  str(zen_deg) + '_' + str(azi_deg) + '.tif'
    array_to_geotiff(DS_max_array, DS_max_path, fin)
    DS_min_path = 'DS_min' + str(cenlat_deg) + '_' +  str(zen_deg) + '_' + str(azi_deg) + '.tif'
    array_to_geotiff(DS_min_array, DS_min_path, fin)
    DS_u0val_path = 'DS_u0val' + str(cenlat_deg) + '_' +  str(zen_deg) + '_' + str(azi_deg) + '.tif'
    array_to_geotiff(DS_u0val_array, DS_u0val_path, fin)
    print("Created DS geotiffs")

    ips_max_path = 'ips_max' + str(cenlat_deg) + '_' +  str(zen_deg) + '_' + str(azi_deg) + '.tif'
    array_to_geotiff(ips_max_array, ips_max_path, fin)
    ips_min_path = 'ips_min' + str(cenlat_deg) + '_' +  str(zen_deg) + '_' + str(azi_deg) + '.tif'
    array_to_geotiff(ips_min_array, ips_min_path, fin)
    ips_u0val_path = 'ips_u0val' + str(cenlat_deg) + '_' +  str(zen_deg) + '_' + str(azi_deg) + '.tif'
    array_to_geotiff(ips_u0val_array, ips_u0val_path, fin)
    print("Created ips geotiffs")

    return PropSumArray, time_kern

def gauss_earth_curvature_radius(center_lat):
    """Calculate Gaussian Earth radius of curvature as a function of latitude.

    center_lat (float): regional latitude in radians
    """
    # Earth ellipse semi-major orbit, REF Wikipedia (https://en.wikipedia.org/wiki/Earth_radius)
    R_equator = 6378.1370 # km (a)

    # Earth ellipse semi-minor axis, REF Wikipedia (https://en.wikipedia.org/wiki/Earth_radius)
    R_polar = 6356.7523142 # km (b)

    # Gaussian Earth radius of curvature, REF Wikipedia (https://en.wikipedia.org/wiki/Earth_radius)
    R_curve = ((R_equator**2.0)*R_polar)/((R_equator*cos(center_lat))**2.0 + (R_polar*sin(center_lat))**2.0)

    # NOT USED PRESENTLY: Directional Earth radius of curvature, REF Wikipedia (https://en.wikipedia.org/wiki/Earth_radius)
    # R_T = (R_polar*R_equator**2)/((R_equator*cos((cent_lat_rad+rel_lat_rad/2)))**2+(R_polar*sin((cent_lat_rad+rel_lat_rad/2)))**2)
    return R_curve

def create_latlon_arrays(R_curve, center_lat, pix_rad):
    """Initial array sizing and create latitude/longitude arrays.

    R_curve (float): radius of earth at latitude
    center_lat (float): regional latitude in radians
    pix_rad (float): pixel size in radians
    """
    p_h = R_curve*pix_rad
    p_w = cos(center_lat)*R_curve*pix_rad

    # Calculate dimensions for kernel
    kernel_cols = int(round(400/p_w))
    kernel_rows = int(round(400/p_h))
    if kernel_cols%2 == 0:
        kernel_cols += 1
    if kernel_rows%2 == 0:
        kernel_rows += 1
    logger.info("kernel width in columns, untrimmed: {}".format(kernel_cols))
    logger.info("kernel height in rows, untrimmed: {}".format(kernel_rows))

    # Create vectors of column and row counts
    col_count = array(range(kernel_cols))
    row_count = array(range(kernel_rows))

    # Find center column, center row
    center_col = int((kernel_cols - 1)/2)
    center_row = int((kernel_rows - 1)/2)

    # Create vecotrs of relative longitude and latitudes (in radians)
    rel_long_vec = array(pix_rad*(col_count - center_col))
    rel_lat_vec = -pix_rad*(row_count - center_row)
    rel_long = tile(rel_long_vec,(kernel_rows,1))
    rel_lat = transpose(tile(rel_lat_vec,(kernel_cols,1)))
    src_lat = rel_lat + center_lat

    return src_lat, rel_long, center_row, center_col

def create_beta(src_lat, rel_long, Chi, azi_view, cen_row, cen_col):
    """Creates array of beta angles (angle between source and observer)

    src_lat (array): array of latitudes
    rel_long (array): array of relative longitudes around 0
    Chi (array): Earth angle from source site
    azi_view (float): viewing azimuth angle in radians
    cen_row (int): row number of center row in arrays
    cen_col (int): column number of center column in arrays
    """
    # beta array, where beta is horizontal angle from direct line (observer to source, OC)
    # to scatter line (observer to scatter, OQ), REF 2, Fig. 6, p. 648

    # formula for horizontal angle on sphere from:
    # http://www.codeguru.com/cpp/cpp/algorithms/article.php/c5115/Geographic-Distance-and-Azimuth-Calculations.htm

    # base trig function
    seterr(invalid='ignore')  # arcsin(NaN) gives invalid warning so ignore it
    beta_arr_raw = arcsin(sin(pi/2.0-src_lat)*sin(rel_long)/sin(Chi))

    # quadrant corrections
    # lower half (left and right)
    beta_arr_raw[cen_row:, :] = pi - beta_arr_raw[cen_row:, :]
    # upper left
    beta_arr_raw[0:cen_row, 0:cen_col] += 2*pi
    # upper right stays the same as original calculation

    return beta_arr_raw - azi_view

def fsum_single(R_T, Chi, u0, l_OC, theta, beta_farg, zen_farg, ubrk_farg, K_am_arg, del_u_farg = .2, lhr_viirs = 1.5):
    """Calculate total light propagation at an observation site along a viewing
    angle vector.
    """
    # initializing dictionary of return values for each pixel of the kernels and debugging outputs for variables within this function
    return_dict = {
        #initializing return variable for kernel pixel value
        "total_sum": 99999,
        # DEBUGGING: initializing return variable for debugging outputs for variables within this function
        "Su_max": 99999,
        "Su_min": 99999,
        "Su_u0_val": 99999,
        "Sd_max": 99999,
        "Sd_min": 99999,
        "Sd_u0_val": 99999,
        "DS_max": 99999,
        "DS_min": 99999,
        "DS_u0_val": 99999,
        "ips_max": 99999,
        "ips_min": 99999,
        "ips_u0_val": 99999
    }
    # if l_OC is null, we can't make this calculation, so return indicator value for null pixels for use in fsum_2d
    if isnan(l_OC):
        return_dict["total_sum"] = 88888
        return return_dict
        #print("L_OC ERROR: R_T:" + str(R_T) + "Chi: " + str(Chi))

    # Scattering distance increment for finite integration over u (km)
    del_u0 = del_u_farg

    # zenith angle at observation site
    zen = zen_farg # radians

    # horizontal angle at observation site from line-of-sight (OC) to line of incoming scattered light (OQ)
    beta = beta_farg #radians

    # length of u, scattering distance, at which the integration increment is relaxed
    ubreak = ubrk_farg #km

    # K - Parameter for relative importance of aerosols to molecules, REF 1, p. 10
    K_am = K_am_arg # unit-less ratio

    #*Constants*
    # N_m,0 - Molecular density at sea level, REF 2, p. 645
    N_m0 = 2.55e19 # cm^-3

    # c - Inverse scale altitude, REF. 2, p. 645
    c_isa = 0.104 # km^-1

    # sigma_m - Integrated Rayleigh scattering cross-section in the V-band, REF 2, p. 646
    sig_m = 1.136e-26 # cm^-2*sr^-1

    # Cinzano shape parameters, REF 2, p. 646
    a1_sp = 0.46
    a2_sp = 0.54
    a3_sp = 0.5 # assumed. No shape factor given in Falchi or Cinzano

    # Falchi best fit parameters for normalzied emssion
    # Weights for angular distributions. REF 1, p. 21
    W_a = 1.9e-3
    W_b = 5.2e-4
    W_c = 7.6e-5
    # light loss per night hour after midnight (negative hours before midnight)
    d_ll = -0.045
    # *End Constants*

    # containers for variables that update in loop
    # DEBUGGING: put in a +0.2 since u_0 seem to be causing the North/South spikes in the kernel values
    u_OQ = u0 + 0.2 # km
    total_sum = 0
    df_prop = 1

    # Total Propagation stable to 3 significant figures
    stability_limit = 0.001

    # loop counter
    lc = 1

    # Upper limit of integration
    # DEBUGGING: u_lim was originally 30.0 but we are running with 250.0 to test if the longer path of
    # integration will help with the longer line of sight towards the horizon which could be
    # influenced by other light sources and atmospheric scattering
    u_lim = 250.0

    # this is where the double summation/integration starts REF 2, pg. 645 equation 9
    while u_OQ < u_lim: # df_prop > stability_limit*total_sum:
        # START OF "u" LOOP

        if u_OQ < ubreak:
            del_u = del_u0
        else:
            del_u = del_u0*10.0

        # s, Distance from source to scattering CQ, REF 2, Appendix A (A1), p. 656
        # equation is wrong in Ref 2 (Cinzano). Changed Chi to theta according to Ref 3, p. 308, Equation 7 (Garstang)
        s_CQ = sqrt((u_OQ - l_OC)**2.0 + 4.0*u_OQ*l_OC*sin(theta/2.0)**2.0) # km

        # h, Height of scattering above Earth reference surface Q, REF 2, Appendix A (A1), p. 656
        h_Q = R_T*(sqrt(1.0 + (u_OQ**2.0 + 2.0*u_OQ*R_T*cos(zen))/R_T**2.0) - 1.0) # km

        # phi, Elevation angle of emission over line of sight from C to O (QCO), REF 2, Appendix A (A1), p. 656
        phi = arccos((l_OC**2.0 + s_CQ**2.0 - u_OQ**2.0)/(2.0*l_OC*s_CQ)) # radians
        phi_deg = phi*180.0/pi # degrees

        # q3, Intermediate quantity, REF 2, Appendix A (A1), p. 656
        q3 = u_OQ*cos(zen)*cos(Chi) - 2.0*R_T*sin(Chi/2.0)**2.0 # km
        # q2, Intermediate quantity (q2=q3 if zenith angle, z=0), REF 2, Appendix A (A1), p. 656
        q2 = u_OQ*sin(zen)*cos(beta)*sin(Chi) + q3 # km

        # Psi, emission angle from source, REF 2, Appendix A (A1), p. 656
        Psi = arccos(q2/s_CQ) # radians
        Psi_deg = Psi*180.0/pi # degrees
        pi_t = 180.0/pi
        # omega, scattering angle at Q, REF 2, Appendix A (A1), p. 656
        omega = theta + phi # radians
        omega_deg = omega*180.0/pi # degrees

        # a, Scale height of aerosols, REF 2, p. 646
        a_sha = 0.657 + 0.059*K_am # km^-1


        if zen > PsiZ_cond:
            # Alternative p2, For zen near 90 degrees, REF 3 p. 309
            p2 = u_OQ - u_OQ**3.0*(8.0/(27.0*pi))*a_sha/R_T

            # Alternative p1, For zen near 90 degrees, REF 3 p. 309
            p1 = u_OQ - u_OQ**3.0*(8.0/(27.0*pi))*c_isa/R_T

        else:
            # p4, Intermediate quantity u-path, REF 2, Appendix A (A2), p. 656
            p4 = (a_sha**2.0 * u_OQ ** 2.0 * cos(zen) ** 2.0 + 2.0 * a_sha * u_OQ * cos(zen) + 2.0) * exp(-a_sha * u_OQ * cos(zen)) - 2.0

            # p3, Intermediate quantity u-path, REF 2, Appendix A (A2), p. 656
            p3 = (c_isa**2.0 * u_OQ ** 2.0 * cos(zen) ** 2.0 + 2.0 * c_isa * u_OQ * cos(zen) + 2.0) * exp(-c_isa * u_OQ * cos(zen)) - 2.0

            # p2, Intermediate quantity u-path, REF 2, Appendix A (A2), p. 656
            # p2 = a_sha**-1.0*cos(zen*(1.0 - exp(-a_sha*u_OQ*cos(zen)) + ((16.0*p4*tan(zen)**2.0)/(9.0*pi*2.0*a_sha*R_T))))**-1.0
            p2 = a_sha**(-1.0)*cos(zen)**(-1.0)*(1.0 - exp(-a_sha*u_OQ*cos(zen)) + ((16.0*p4*tan(zen)**2.0)/(9.0*pi*2.0*a_sha*R_T)))

            # p1, Intermediate quantity u-path, REF 2, Appendix A (A2), p. 656
            # p1 = c_isa**-1.0*cos(zen*(1.0 - exp(-c_isa*u_OQ*cos(zen)) + ((16.0*p3*tan(zen)**2.0)/(9.0*pi*2.0*c_isa*R_T))))**-1.0
            p1 = c_isa**(-1.0)*cos(zen)**(-1.0)*(1.0 - exp(-c_isa*u_OQ*cos(zen)) + ((16.0*p3*tan(zen)**2.0)/(9.0*pi*2.0*c_isa*R_T)))

        # ksi1, Extinction of light along u-path from scatter at Q to observation at O, REF 2, Appendix A (A2), p. 656
        ksi1 = exp(-N_m0*sig_m*(p1 + 11.778*K_am*p2))

        if Psi > PsiZ_cond:
            # Alternative f2, For zen near 90 degrees, REF 3 p. 309
            f2 = s_CQ - s_CQ**3.0*(8.0/(27.0*pi))*a_sha/R_T

            # Alternative f1, For zen near 90 degrees, REF 3 p. 309
            f1 = s_CQ - s_CQ**3.0*(8.0/(27.0*pi))*c_isa/R_T

        else:
            # f4, Intermediate quantity s-path, REF 2, Appendix A (A2), p. 657
            f4 = (a_sha**2.0*s_CQ**2.0*cos(Psi)**2.0 + 2.0*a_sha*s_CQ*cos(Psi) + 2.0)*exp(-a_sha*s_CQ*cos(Psi)) - 2.0

            # f3, Intermediate quantity s-path, REF 2, Appendix A (A2), p. 657
            f3 = (c_isa**2.0*s_CQ**2.0*cos(Psi)**2.0 + 2.0*c_isa*s_CQ*cos(Psi) + 2.0)*exp(-c_isa*s_CQ*cos(Psi)) - 2.0

            # f2, Intermediate quantity s-path, REF 2, Appendix A (A2), p. 657
            # f2 = a_sha**-1.0*cos(Psi*(1.0 - exp(-a_sha*s_CQ*cos(Psi)) + ((16.0*f4*tan(Psi)**2.0)/(9.0*pi*2.0*a_sha*R_T))))**-1.0
            f2 = a_sha**(-1.0)*cos(Psi)**(-1.0)*(1.0 - exp(-a_sha*s_CQ*cos(Psi)) + ((16.0*f4*tan(Psi)**2.0)/(9.0*pi*2.0*a_sha*R_T)))

            # f1, Intermediate quantity s-path, REF 2, Appendix A (A2), p. 657
            # f1 = c_isa**-1.0*cos(Psi*(1.0 - exp(-c_isa*s_CQ*cos(Psi)) + ((16.0*f3*tan(Psi)**2.0)/(9.0*pi*2.0*c_isa*R_T))))**-1.0
            f1 = c_isa**(-1.0)*cos(Psi)**(-1.0)*(1.0 - exp(-c_isa*s_CQ*cos(Psi)) + ((16.0*f3*tan(Psi)**2.0)/(9.0*pi*2.0*c_isa*R_T)))

        # ksi2, Extinction of light along s-path from emission at C to scatter at Q, REF 2, Appendix A (A2), p. 656
        ksi2 = exp(-N_m0*sig_m*(f1 + 11.778*K_am*f2))

        # I(Psi), Normalized emission function, MODIFIED FROM REF 1, p. 13 (leaving out natural sky brightness)
        # There is no indication in References of a shape parameter for the W_c term of the equation
        I_ne = 1.0/(2.0*pi)*(W_a*2.0*a1_sp*cos(Psi) + W_b*0.554*a2_sp*Psi**4.0 + W_c*a3_sp*sin(3.0*Psi))*(1.0+lhr_viirs*d_ll)
        # i(Psi,s), Illuminance per unit flux, REF 2, Eq. 6, p. 644
        i_ps = I_ne*ksi2/s_CQ**2.0

        # N_m(h), Number density of gaseous component of atmosphere as function of altitude, h, REF 2, Eq. 10, p. 645
        N_m = N_m0*exp(-c_isa*h_Q)

        # N_a*sigma_a, Total integrated scattering cross-section, REF 2, Eq. 12, p. 645 **REARRANGED**
        Na_x_siga = K_am*N_m*sig_m*11.11

        # f_m, Angular scattering function for molecular Rayleigh scattering, REF 2, Eq. 13, p. 646
        f_m = 3.0*(1.0 + cos(omega)**2.0)/(16.0*pi)

        # f_a, Angular scattering function for aerosol Mie scattering, REF 2, Eq. 14, p. 646
        if 0.0 <= omega_deg <= 10.0:
            f_a = 7.5*exp(-0.1249*omega_deg**2.0/(1.0 + 0.04996*omega_deg**2.0))
        elif 10.0 < omega_deg <= 124.0:
            f_a = 1.88*exp(-0.07226*omega_deg + 0.0002406*omega_deg**2.0)
        elif 124.0 < omega_deg <= 180.0:
            f_a = 0.025 + 0.015*sin((2.25*omega_deg - 369.0)*(pi/180.0))

        # S_d, Luminous flux per unit solid angle per unit upward flux (directly from source), REF 2, Eq. 5, p. 644
        S_d = (N_m*sig_m*f_m + Na_x_siga*f_a)*i_ps

        # D_S, Double scattering correction factor, REF 2, Eq. 20, p. 647
        D_S = 1 + N_m0*sig_m*(11.11*K_am*f2 + (f1/3.0))

        # S_u, Total illumance as a function of u, REF 2, Eq. 8, p. 645
        S_u = S_d*D_S

        df_prop = S_u*ksi1*del_u
        # integrand of propogation function, REF 2, Eq. 3, p. 644
        total_sum += df_prop
        u_OQ += del_u

        # DEBUGGING: assign debugging values
        # for each value we are debugging, we locate the maximum, minimum, and
        # u_0 starting value. We check each debugging variable to check if it's at
        # the initializing value of 99999, if so assign the current variable value,
        # otherwise compare the variable to the existing max or min to see if then
        # max or min needs to be updated
        # DEBUGGING: S_u debugging returns
        if abs(return_dict["Su_max"]- 99999) < 0.1:
            return_dict["Su_max"] = S_u
            return_dict["Su_u0_val"] = S_u
        else:
            if S_u > (return_dict["Su_max"]):
                return_dict["Su_max"] = S_u

        if abs(return_dict["Su_min"] - 99999) < 0.1:
            return_dict["Su_min"] = S_u
        else:
            if S_u < (return_dict["Su_min"]):
                return_dict["Su_min"] = S_u

        # DEBUGGING: Sd debugging returns
        if abs(return_dict["Sd_max"]- 99999) < 0.1:
            return_dict["Sd_max"] = S_d
            return_dict["Sd_u0_val"] = S_d
        else:
            if S_d > (return_dict["Sd_max"]):
                return_dict["Sd_max"] = S_d

        if abs(return_dict["Sd_min"] - 99999) < 0.1:
            return_dict["Sd_min"] = S_d
        else:
            if S_d < (return_dict["Sd_min"]):
                return_dict["Sd_min"] = S_d

        # DEBUGGING: DS debugging returns
        if abs(return_dict["DS_max"]- 99999) < 0.1:
            return_dict["DS_max"] = D_S
            return_dict["DS_u0_val"] = D_S
        else:
            if D_S > (return_dict["DS_max"]):
                return_dict["DS_max"] = D_S

        if abs(return_dict["DS_min"] - 99999) < 0.1:
            return_dict["DS_min"] = D_S
        else:
            if D_S < (return_dict["DS_min"]):
                return_dict["DS_min"] = D_S

        # DEBUGGING:i_ps debugging returns
        if abs(return_dict["ips_max"]- 99999) < 0.1:
            return_dict["ips_max"] = i_ps
            return_dict["ips_u0_val"] = i_ps
        else:
            if i_ps > (return_dict["ips_max"]):
                return_dict["ips_max"] = i_ps

        if abs(return_dict["ips_min"] - 99999) < 0.1:
            return_dict["ips_min"] = i_ps
        else:
            if i_ps < (return_dict["ips_min"]):
                return_dict["ips_min"] = i_ps

    if total_sum > 1:
        print('break!')
    # assign the pixel value for this pixel of the kernel, note that this Value
    # should be very small, generally in the range of 10^(-16) to 10^(-13)
    return_dict["total_sum"] = total_sum

    return return_dict

def array_to_geotiff(array, outfilename, referenceVIIRS, new_trans = None):
    """Save numpy array as a geotiff.

    array (array): numpy array to save
    outfilename (str): output file name
    referenceVIIRS (str): path to VIIRS DNB tif to extract tif info
    new_trans (list): GDAL geo transform tuple (in list form) (default None)
    """
    imdata = gdal.Open(referenceVIIRS)

    # Save out to a GeoTiff
    arr = array
    # First, gather some information from VIIRS image
    [cols,rows] = arr.shape
    # set transform to VIIRS transfomr if there is no new transform
    trans = imdata.GetGeoTransform() if not new_trans else new_trans
    proj = imdata.GetProjection()
    nodatav = -1000000000
    outfile = outfilename
    logger.info('Saving geotiff {}'.format(outfile))

    # Create the georeffed file, using the information from the VIIRS image
    outdriver = gdal.GetDriverByName("GTiff")
    outdata = outdriver.Create(str(outfile), rows, cols, 1, gdal.GDT_Float32)

    # Write data to the file, which is the kernel array in this example
    outdata.GetRasterBand(1).WriteArray(arr)

    # Set a no data value
    outdata.GetRasterBand(1).SetNoDataValue(nodatav)

    # Georeference the image
    outdata.SetGeoTransform(trans)

    # Write projection information
    outdata.SetProjection(proj)

if __name__ == "__main__":
    logging.basicConfig(level=logging.DEBUG, format='%(levelname)s: %(message)s')
    main()
