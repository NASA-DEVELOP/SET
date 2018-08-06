import argparse
import sys
import logging
logger = logging.getLogger()
from multiprocessing import Pool

from numpy import loadtxt

import skyglow.skyglow as skyglow
import skyglow.darksky as darksky


def main():
    logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
    parser = argparse.ArgumentParser(description='CLI for calculating skyglow')

    # generic args
    parser.add_argument('action', choices=['sgmap_single', 'kernel_lib', 'sgmap_multiple', 'gui', 'hemisphere'])
    parser.add_argument('-l', '--latitude', type=float, help='Site/hemishpere center latitude')
    parser.add_argument('viirs_file', type=str, nargs='?', help='VIIRS image file path')
    parser.add_argument('--debug', help='increase output verbosity', action='store_true')

    # sgmap_single specific args
    parser.add_argument('-kam', '--clarity', type=float, default=1, help='Atmospheric clarity ratio (default: 1)')
    parser.add_argument('-z', '--zenith', type=float, default=0, help='Zenith angle (degrees) (default: 0)')
    parser.add_argument('-a', '--azimuth', type=float, default=0, help='Azimuth angle (degrees) (default: 0)')
    parser.add_argument('-k', '--kernel_file', type=str, default='', help='sgmap_single kernel file path')

    # kernel_lib specific args
    parser.add_argument('-ang', '--angles_csv', type=str, default='default.csv', help='kernel_lib angles CSV path (default: default.csv)')
    parser.add_argument('-hem', '--hemisphere', action='store_true', help='Generate kernels for hemispherical visualization')
    parser.add_argument('-s', '--sync', action='store_true', help='Generate kernels concurrently or not')

    # sgmap_multiple specific args
    parser.add_argument('-klib', '--kernel_folder', type=str, help='sgmap_multiple kernel folder path')
    parser.add_argument('-out', '--output_folder', type=str, help='sgmap_multiple output folder path')

    # hemisphere specific args
    parser.add_argument('-sglib', '--skyglow_folder', type=str, help='Hemisphere skyglow map folder path')
    parser.add_argument('-lon', '--longitude', type=float, help='Hemisphere center longitude')


    args = parser.parse_args()
    if args.debug:
        logger.setLevel(logging.DEBUG)
    if args.action == 'gui':
        skyglow.main()
    elif args.action == 'hemisphere':
        if not args.latitude:
          raise ValueError('Latitude is required for hemisphere')
        if not args.longitude:
          raise ValueError('Longitude is required for hemisphere')
        if not args.skyglow_folder:
          raise ValueError('Skyglow folder path is required for hemisphere')
        darksky.generate_hem(args.latitude, args.longitude, args.skyglow_folder)
    else:
        if not args.viirs_file:
            raise ValueError('VIIRS file path is required for sgmap_single, kernel_lib, sgmap_multiple actions')

        if args.action == 'sgmap_single':
            if not args.latitude:
                raise ValueError('Latitude is required for sgmap_single')
            if args.kernel_file:
                print('Using kernel file {} for sgmap_single'.format(args.kernel_file))
                krn_file = args.kernel_file
            darksky.sgmapper(args.latitude, args.clarity, args.zenith, args.azimuth, args.viirs_file, args.kernel_file)
        elif args.action == 'kernel_lib':
            if not args.latitude:
                raise ValueError('Latitude is required for kernel_lib')
            if not args.sync:
                static_args = [(args.latitude, args.clarity, args.viirs_file, args.hemisphere)]
                with open(args.angles_csv, "rb") as f:
                    angle_list = loadtxt(f, delimiter=",",skiprows=1)
                # cartesian product of arguments
                args_product = [(x[0], x[1], y[0], y[1], x[2], x[3]) for x in static_args for y in angle_list]
                p = Pool()
                try:
                    p.map_async(darksky.krn_unpack, args_product).get(999999)
                except KeyboardInterrupt:
                    p.terminate()
            else:
                darksky.krnlibber(args.latitude, args.clarity, args.angles_csv, args.viirs_file, args.hemisphere)
        elif args.action == 'sgmap_multiple':
            if not args.kernel_folder:
                raise ValueError('Kernel folder path is required for sgmap_multiple')
            if not args.output_folder:
                raise ValueError('Output folder path is required for sgmap_multiple')
            darksky.multisgmapper(args.viirs_file, args.kernel_folder, args.output_folder)

if __name__ == '__main__':
    main()
