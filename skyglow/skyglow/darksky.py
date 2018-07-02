# REFERENCES
# (1) Falchi, F., P. Cinzano, D. Duriscoe, C.C.M. Kyba, C.D. Elvidge, K. Baugh,
#   B.A. Portnov, N.A. Rybnikova and R. Furgoni, 2016. The new workd atlas of
#   artificial night sky brightness. Sci. Adv. 2.
# (2) Cinzano, P., F. Falchi, C.D. Elvidge and  K.E. Baugh, 2000. The
#   artificial night sky brightness mapped from DMSP satellite Operational
#   Linescan System measurements. Mon. Not. R. Astron. Soc. 318.
# (3) Garstang, R.H., 1989. Night-sky brightness at observatories and sites.
#   Pub. Astron. Soc. Pac. 101.

from __future__ import division
from __future__ import print_function

from numpy import *
# https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.fftconvolve.html
# scipy.signal.fftconvolve
# from scipy import ndimage
from matplotlib import pyplot as plt
import matplotlib.colors as colors
from osgeo import gdal

import constants

import os.path, ntpath
import time
import json
import logging
import sys
import glob
import warnings
from multiprocessing import Pool
from builtins import range
warnings.filterwarnings("error")
logger = logging.getLogger()


default_kernel_folder = os.path.join(os.getcwd(), "kernel_lib")
default_output_folder = os.path.join(os.getcwd(), "skyglow_out")

# ISSUES
# check angle fixes ~sgmapper and other
# check beta array fixes ~line 298
# continue abstracting prop2d is sgmapper ~line 50
# still need to get rid of global variables
# still consider making external variable/function constants file
# consider whether name "darkskypy" is okay
# darkskypy mixes spaces and tabs for indentation

PsiZ_cond = 89.5*(pi/180)
ubr_arg = 10.0

logger.info(sys.version)

def main():
    """Main entry point if darksky is run from the command line. Calls
       appropriate action."""

    if len(sys.argv) < 2:
        raise ValueError("Usage: python {} <action>. Action can be 'sgmap_single', 'kernel_lib', 'sgmap_multiple'".format(sys.argv[0]))

    project_root = os.path.join(os.getcwd(), '/'.join(sys.argv[0].split("/")[:-1]))

    action = sys.argv[1]

    # Estimate a single map of artificial sky brightness with or without
    # an existing 2d propagation function or "kernel"
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

    # Estimate a group or "library" of 2d propagation functions or "kernels"
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
            sync = bool(sys.argv[5])
        except:
            sync = False

        if not sync:
            static_args = [(la, ka, fi)]
            with open(csv_path, "rb") as f:
                angle_list = loadtxt(f, delimiter=",",skiprows=1)
            # cartesian product of arguments
            args_product = [(x[0], x[1], x[2], y[0], y[1]) for x in static_args for y in angle_list]
            p = Pool()
            try:
                p.map_async(krn_unpack, args_product).get(9999999)
            except KeyboardInterrupt:
                p.close()
        else:
            krnlibber(la, ka, csv_path, fi)

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

    else:
        raise ValueError("No action/incorrect action given. Choose ""sgmap_single"", ""kernel_lib"", or ""sgmap_multiple""")

def sgmapper(centerlat_arg, k_am_arg, zen_arg, azi_arg, filein, prop2filein=""):
    kerneltiffpath = prop2filein
    # added code to allow for auto opening of kernels based on if they are saved in the current directory or the kernel_lib directory
    if kerneltiffpath == "":
        # If there is no input from skyglow.py for in input kernel (such as when the program is run through command line) then the default kernel name is used
        kerneltiffpath = 'kernel_' + str(centerlat_arg) + '_' +  str(k_am_arg) + '_' + str(zen_arg) + '_' + str(azi_arg) + '.tif'
    if (os.path.isfile(kerneltiffpath) is False) and (os.path.isfile(os.path.join("kernel_lib", kerneltiffpath)) is False):
        # If there is no 2d propagation function (kernel), ...
        # then estimate the 2d propagation function

        # first convert angles in degrees to radians
        azi_rad = azi_arg*(pi/180)
        zen_rad = zen_arg*(pi/180)
        if abs(centerlat_arg) > 90:
            raise ValueError("Latitude must be between 0 and 90 degrees")
        lat_rad = centerlat_arg*(pi/180)
        propkernel, totaltime = fsum_2d(lat_rad, k_am_arg, zen_rad, azi_rad, filein)
        kerneltiffpath = 'kernel_' + str(centerlat_arg) + '_' +  str(k_am_arg) + '_' + str(zen_arg) + '_' + str(azi_arg) + '.tif'
        array_to_geotiff(propkernel, kerneltiffpath, filein)

        logger.info("time for prop function ubreak 10: %s", totaltime)
    else:
        # If there is a 2d propagation function (kernel) in the current direcotry and not in the kernel_lib subdirectory, ...
        # then read it into an array
        if (os.path.isfile(kerneltiffpath) is True) and (os.path.isfile(os.path.join("kernel_lib", kerneltiffpath)) is False):
            kerneldata = gdal.Open(kerneltiffpath)
            propkernel = kerneldata.ReadAsArray()
        # If there is a 2d propagation function (kernel) in the kernel_lib sub directory but not the current directory, ...
        # then read it into an array
        elif (os.path.isfile(os.path.join("kernel_lib", kerneltiffpath)) == True) and (os.path.isfile(kerneltiffpath) == False):
            kerneldata = gdal.Open(os.path.join("kernel_lib", kerneltiffpath))
            propkernel = kerneldata.ReadAsArray()
        # If there is a 2d propagation function (kernel) in the current directory and subdirectory, ...
        # then read the one in the current directory into an array
        elif (os.path.isfile(kerneltiffpath) == True) and (os.path.isfile(os.path.join("kernel_lib",kerneltiffpath)) == True):
            kerneldata = gdal.Open(kerneltiffpath)
            propkernel = kerneldata.ReadAsArray()

    # read in VIIRS DNB image
    viirsraster = gdal.Open(filein)
    imagearr = viirsraster.ReadAsArray()

    # ADDED DEBUG CODE #
    logger.debug('imagearr shape : {}'.format(imagearr.shape))
    logger.debug('imagearr shape rows : {}'.format(imagearr.shape[0]))
    logger.debug('imagearr shape columns : {}'.format(imagearr.shape[1]))
    logger.debug('propkernel shape : {}'.format(propkernel.shape))
    logger.debug('propkernel shape rows : {}'.format(propkernel.shape[0]))
    logger.debug('propkernel shape columns : {}'.format(propkernel.shape[1]))

    # Produce sky glow raster
    skyglowarr = convolve_viirs_to_skyglow(imagearr, propkernel)
    skyglowpath = (ntpath.basename(filein)[:-4] + '_' + str(centerlat_arg) + '_' + str(ubr_arg) + '_'
                   + str(zen_arg) + '_' + str(azi_arg) + 'convolved.tif')
    array_to_geotiff(skyglowarr, skyglowpath, filein)
    logger.info("===============\n***Finished!***\n===============\nSkyglow Map saved as:\n" + skyglowpath)
    constants.ding()

def krn_unpack(args):
    print(args)
    return generate_krn(*args)

def generate_krn(centerlat_arg, k_am_arg, zen_arg, azi_arg, filein):
    lat_rad = centerlat_arg*(pi/180)
    zen_rad = zen_arg*(pi/180)
    azi_rad = azi_arg*(pi/180)
    propkernel, totaltime = fsum_2d(lat_rad, k_am_arg, zen_rad, azi_rad, filein)
    kerneltiffpath = 'kernel_' + str(centerlat_arg) + '_' +  str(k_am_arg) + '_' + str(zen_arg) + '_' + str(azi_arg) + '.tif'
    array_to_geotiff(propkernel, kerneltiffpath, filein)
    logger.info("time for prop function ({}, {}) ubreak 10: %s".format(zen_arg, azi_arg), totaltime)

def krnlibber(centerlat_arg, k_am_arg, angles_file, filein):
    lat_rad = centerlat_arg*(pi/180)
    with open(angles_file, "rb") as f:
        angle_list = loadtxt(f, delimiter=",",skiprows=1)
    for angle_set in angle_list:
        zenith = angle_set[0]
        azimuth = angle_set[1]
        zen_rad = zenith*(pi/180)
        azi_rad = azimuth*(pi/180)

        propkernel, totaltime = fsum_2d(lat_rad, k_am_arg, zen_rad, azi_rad, filein)
        kerneltiffpath = 'kernel_' + str(centerlat_arg) + '_' +  str(k_am_arg) + '_' + str(zenith) + '_' + str(azimuth) + '.tif'
        array_to_geotiff(propkernel, kerneltiffpath, filein)
        logger.info("time for prop function ({}, {}) ubreak 10: %s".format(zenith, azimuth), totaltime)

def multisgmapper(filein, krnfolder, outfolder):
    # read in VIIRS DNB image
    viirsraster = gdal.Open(filein)
    imgarr = viirsraster.ReadAsArray()

    # Get kernel library file list
    krnsrch = os.path.join(krnfolder,'*.tif')
    krnlist = glob.glob(krnsrch)

    ln = 0

    for krn in krnlist:
        # If the loop has already ran once then rescale imgarr back down to prevent multiple scaling factors being applied
        if ln > 0:
            imgarr *= 10**-9
        ln = 1

        # If there is a 2d propagation function (kernel), ...
        # then read it into an array
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
        print(sgfile)
        # Create skyglow array
        sgarr = convolve_viirs_to_skyglow(imgarr, krnarr)
        # Write skyglow to a geotiff
        array_to_geotiff(sgarr, sgfile, filein)
        # Tests for an input azimuth angle of between 0 and 180 degrees and does the complimenting kernel as well
        cond_angle = sgtags[4]
        cond = float(cond_angle[:-4])
        if (cond < 180.0) and not (cond == 0.0):
            # Must undo the scaling factor in convolve_viirs_to_skyglow
            rescale_factor = 10**-9
            img_rescale = rescale_factor * imgarr
            # Flips the current kernel left to right
            krnarr_flip = fliplr(krnarr)
            # Creates the corresponding flipped angle
            flip_azi = 360.0 - cond
            # Returns the proper file extention onto the tag in order to read as a file
            sgtags[4] = str(flip_azi) + '.tif'
            sgbase_flip = sep.join(sgtags)
            sgfile_flip = os.path.join(outfolder, sgbase_flip)
            print(sgfile_flip)
            logger.info("working on " + sgbase_flip)
            sgarr_flip = convolve_viirs_to_skyglow(img_rescale, krnarr_flip)
            array_to_geotiff(sgarr_flip, sgfile_flip, filein)

# Function convolves VIIRS DNB with 2d propagation function
def convolve_viirs_to_skyglow(dnbarr, proparr):
    # Convert to float 32 for Fourier, scale, and round
    # falchi assumed natural sky brightness to be 174 micro cd/m^2 = 2.547e-11 watt/cm^2/steradian (at 555nm)
    # Not sure if this is correct scaling factor, I assume that this makes the output prop image in units of cd/m^2
    viirs_scaling_factor = 10**9
    dnbarr *= viirs_scaling_factor
    proparr = float32(nan_to_num(proparr))

    # ADDED DEBUG CODE
    logger.debug('dnbarr.shape[0] : {}'.format(dnbarr.shape[0]))
    logger.debug('dnbarr.shape[1] : {}'.format(dnbarr.shape[1]))
    logger.debug('proparr.shape[0] : {}'.format(proparr.shape[0]))
    logger.debug('proparr.shape[1] : {}'.format(proparr.shape[1]))

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

    # ADDED DEBUG CODE
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

    # Fourier Transform Method #

    np_dft_prop_im = fft.fft2(padded_prop)
    np_dft_kernel_shift = fft.fftshift(np_dft_prop_im)
    np_magnitude_spectrum = 20*log(abs(np_dft_kernel_shift))

    # ADDED DEBUG CODE
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

# Function that creates 2d propagation function
def fsum_2d(cenlat, k_am, zen, azi, fin):
    # Input Variables
    logger.info('**INPUTS**')

    # 2d propagation function pixel size in degrees and radians
    p_deg = .0041666667
    p_rad = p_deg*pi/180.0

    # central latitude for estimation of the 2d propagation function
    logger.info('Central latitude (deg): {}'.format(cenlat*(180/pi)))

    # central latitude for estimation of the 2d propagation function
    logger.info('Atmospheric clarity indicator: {}'.format(k_am))

    # z, Zenith angle site, REF 2, Fig. 6, p. 648
    logger.info('z, Site zenith (deg): {}'.format(zen*(180/pi)))

    # Beta, Azimuth angle site, REF 2, Fig 6 p. 648
    logger.info('Beta, Site azimuth as referenced from North (deg): {}'.format(azi*(180/pi)))

    # ubr, Length of u for relaxing integration increment
    ubr = ubr_arg
    logger.info('ubr, Length of u for relaxing integration increment (km): {}'.format(ubr))

    # Gaussian Earth radius of curvature, REF Wikipedia (https://en.wikipedia.org/wiki/Earth_radius)
    R_T = gauss_earth_curvature_radius(cenlat)
    logger.info('R_T, Radius of curvature of the Earth (km): {}'.format(R_T))

    # Create latitude/longitude arrays
    source_lat, rltv_long, cent_row, cent_col = create_latlon_arrays(R_T, cenlat, p_rad)

    # Distance from source (C) to observation site (O) along ellipsoid surface, REF 2, Fig. 6, p. 648
    # using haversine formula
    D_OC = 2.0*R_T*arcsin(sqrt(sin((source_lat - cenlat)/2.0)**2.0 + cos(cenlat)*cos(source_lat)*sin(rltv_long/2.0)**2.0))

    # Assignment of NaNs or null values outside of 200 km
    D_OC[D_OC > 201.0] = NaN

    # Check of D_OC shape after assigning NaNs outside of 200 km radius
    # WHAT IS THE POINT OF CHECKING THE SHAPE IF WE DIDN'T CHANGE IT!
    logger.info('kernel heigth in pixels, trimmed: {}'.format(D_OC.shape[0]))
    logger.info('kernel width in pixels, trimmed: {}'.format(D_OC.shape[1]))
    widthcenter = (D_OC.shape[1] + 1)//2
    heigthcenter = (D_OC.shape[0] + 1)//2
    # reassignment of center value, need to use better method
    D_OC[cent_row,cent_col] = .01

    # Earth angle from source to site, REF 3, p. 308\
    Chi = D_OC/R_T

    beta = create_beta(source_lat, rltv_long, Chi, azi, cent_row, cent_col)

    # u0, shortest scattering distance based on curvature of the Earth, REF 2, Eq. 21, p. 647
    u0 = 2.0*R_T*sin(Chi/2.0)**2.0/(sin(zen)*cos(beta)*sin(Chi)+cos(zen)*cos(Chi)) #km
    #array_to_geotiff(u0, "u0asymmetrytest.tif", fin)

    # l, Direct line of sight distance between source and observations site, REF 2, Appendix A (A1), p. 656
    # l_OC and D_OC are similar as expected
    l_OC = sqrt(4.0*R_T**2.0*sin(Chi/2.0)**2.0) # km

    # q1, Intermediate quantity, REF 2, Appendix A (A1), p. 656, **WITH CORRECTION FROM REF 3, eq. 6, p. 308**
    q1 = R_T*(sin(Chi)*sin(zen)*cos(beta) + cos(Chi)*cos(zen) - cos(zen)) # km

    # theta, elevation angle of scatter above source from site (QOC), REF 2, Appendix A (A1), p. 656
    theta = arccos(q1/l_OC) # radians

    # simulates what a 180 degrees azimuth input would look like to the program to make sure that rounding errors dont prevent the if statement
    logger.info('azi, {}'.format(azi))
    azi_of_180 = float(180)*pi/180.0
    if azi in {0.0, azi_of_180}:
        # Get left arrays to cut processing time in half
        Chileft = Chi[0:, 0:widthcenter]
        u0left = u0[0:, 0:widthcenter]
        l_OCleft = l_OC[0:, 0:widthcenter]
        thetaleft = theta[0:, 0:widthcenter]
        betaleft = beta[0:, 0:widthcenter]
        PropSumArrayleft = zeros_like(l_OCleft)
        kerdim = l_OCleft.shape
        ker10per = kerdim[0]//10

        logger.info("Beginning kernel...")
        start = time.time()

        # 2d iteration for integrating from u0 to infinity to create propagation function for each element
        #for p,c,u,l,t in itertools.izip(nditer(PropSumArrayleft, op_flags=['readwrite']),nditer(Chileft, op_flags=['readwrite']),nditer(u0left, op_flags=['readwrite']), nditer(l_OCleft, op_flags=['readwrite']),nditer(thetaleft, op_flags=['readwrite'])):
        #    p[...] = fsum_single(R_T, c, u, l, t, 0.0, zen, ubr, k_am)

        for ii in range(kerdim[0]):
            if mod(ii, ker10per) == 0:
                interm = time.time() - start
                logger.info("Kernel (%d, %d), %d percent complete, %d minutes", zen, azi, 100.0*ii/kerdim[0], interm/60.0)
            for jj in range(kerdim[1]):
                PropSumArrayleft[ii][jj] = fsum_single(R_T, Chileft[ii][jj],
                                                       u0left[ii][jj],
                                                       l_OCleft[ii][jj],
                                                       thetaleft[ii][jj],
                                                       betaleft[ii][jj],
                                                       zen, ubr, k_am)

        time_kern = time.time() - start
        logger.info("Time to produce kernel: %d minutes", time_kern/60.0)

        # Mirror left side of kernel array to the right
        PropSumArrayright = fliplr(PropSumArrayleft[:,1:])

        # Complete 2d propagation function by putting left and right together
        PropSumArray = hstack((PropSumArrayleft, PropSumArrayright))
    else:
        # debug checks for subset arrays
        logger.debug('Problem Index Chi: %s', Chi[523, 470])
        logger.debug('Chi: %s', Chi)
        logger.debug('Problem Index u0: %s', u0[523, 470])
        logger.debug('u0: %s', u0)
        logger.debug('Problem Index l_OC: %s', l_OC[523, 470])
        logger.debug('l_OC: %s', l_OC)
        logger.debug('Problem Index theta: %s', theta[523, 470])
        logger.debug('theta: %s', theta)
        logger.debug('Surrounding Indices beta: %s', beta[522:525, 470])
        logger.debug('Problem Index beta: %s', beta[523, 470])
        logger.debug('beta: %s', beta)

        # initialize for 2d iteration
        PropSumArray = zeros_like(l_OC)
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
                logger.info("Kernel (%d, %d), %d percent complete, %d minutes", zen, azi, 100.0*ii/kerdim[0], interm/60.0)
            for jj in range(kerdim[1]):
                PropSumArray[ii][jj] = fsum_single(R_T, Chi[ii][jj], u0[ii][jj],
                                                   l_OC[ii][jj], theta[ii][jj],
                                                   beta[ii][jj], zen, ubr, k_am)

        time_kern = time.time() - start
        logger.info("Time to produce kernel: %d minutes", time_kern/60.0)

        logger.debug('Surrounding Indices PropSumArray: %s', PropSumArray[522:525, 470])
        logger.debug('Problem Index PropSumArray: %s', PropSumArray[523, 470])
        logger.debug('PropSumArray: %s', PropSumArray)

    return PropSumArray, time_kern

# Function to calculate Gaussian Earth radius of curvature as a function of latitude
def gauss_earth_curvature_radius(center_lat):
    # Earth ellipse semi-major orbit, REF Wikipedia (https://en.wikipedia.org/wiki/Earth_radius)
    R_equator = 6378.1370 # km (a)

    # Earth ellipse semi-minor axis, REF Wikipedia (https://en.wikipedia.org/wiki/Earth_radius)
    R_polar = 6356.7523142 # km (b)

    # Gaussian Earth radius of curvature, REF Wikipedia (https://en.wikipedia.org/wiki/Earth_radius)
    R_curve = ((R_equator**2.0)*R_polar)/((R_equator*cos(center_lat))**2.0 + (R_polar*sin(center_lat))**2.0)

    # NOT USED PRESENTLY: Directional Earth radius of curvature, REF Wikipedia (https://en.wikipedia.org/wiki/Earth_radius)
    # R_T = (R_polar*R_equator**2)/((R_equator*cos((cent_lat_rad+rel_lat_rad/2)))**2+(R_polar*sin((cent_lat_rad+rel_lat_rad/2)))**2)

    return R_curve

# Function does initial array sizing and create latitude/longitude arrays
def create_latlon_arrays(R_curve, center_lat, pix_rad):
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

# Function that creates array of beta angles based on
# source latitude array, relative longitude array, Earth angle, and input viewing azimuth
def create_beta(src_lat, rel_long, Chi, azi_view, cen_row, cen_col):
    # beta array, where beta is horizontal angle from direct line (observer to source, OC)
    # to scatter line (observer to scatter, OQ), REF 2, Fig. 6, p. 648

    # formula for horizontal angle on sphere from:
    # http://www.codeguru.com/cpp/cpp/algorithms/article.php/c5115/Geographic-Distance-and-Azimuth-Calculations.htm

    # base trig function
    seterr(invalid='ignore')  # arcsin(NaN) gives invalid warning so ignore it
    beta_arr_raw = arcsin(sin(pi/2.0-src_lat)*sin(rel_long)/sin(Chi))

    # quadrant corrections
    # lower half (left and right)
    beta_arr_raw[0:cen_row, :] = pi - beta_arr_raw[0:cen_row, :]
    # upper left
    beta_arr_raw[cen_row:, 0:cen_col] += 2*pi
    # upper right stays the same as original calculation

    return beta_arr_raw - azi_view

# Function that takes elements of the arrays of D_OC, Chi, etc. as arguments.
def fsum_single(R_T, Chi, u0, l_OC, theta, beta_farg, zen_farg, ubrk_farg, K_am_arg, del_u_farg = .2, lhr_viirs = 1.5):
    if isnan(l_OC):
        return nan

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
    u_OQ = u0 # km
    total_sum = 0
    df_prop = 1

    # Total Propagation stable to 3 significant figures
    stability_limit = 0.001

    # loop counter
    lc = 1

    while u_OQ < 30.0: # df_prop > stability_limit*total_sum:
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
    if total_sum > 1:
        print('break!')
    return total_sum

def array_to_geotiff(array, outfilename, referenceVIIRS):
    imdata = gdal.Open(referenceVIIRS)

    # Save out to a GeoTiff
    arr = array
    # First of all, gather some information from VIIRS image
    [cols,rows] = arr.shape
    trans = imdata.GetGeoTransform()
    proj = imdata.GetProjection()
    nodatav = 0
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
