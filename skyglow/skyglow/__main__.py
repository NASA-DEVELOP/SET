import argparse
import sys
import logging
logger = logging.getLogger()
from multiprocessing import Pool
from itertools import product

from numpy import loadtxt

import skyglow
import darksky


def main():
    logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
    parser = argparse.ArgumentParser(description='CLI for calculating skyglow')
    parser.add_argument('action', choices=['sgmap_single', 'kernel_lib', 'sgmap_multiple', 'gui'])
    parser.add_argument('-l', '--latitude', type=float, help='Site latitude')
    parser.add_argument('-kam', '--clarity', type=float, default=1, help='Atmospheric clarity ratio (default: 1)')
    parser.add_argument('-z', '--zenith', type=float, default=0, help='Zenith angle (degrees) (default: 0)')
    parser.add_argument('-a', '--azimuth', type=float, default=0, help='Azimuth angle (degrees) (default: 0)')
    parser.add_argument('-k', '--kernel_file', type=str, help='sgmap_single kernel file path')

    parser.add_argument('-ang', '--angles_csv', type=str, default='default.csv', help='kernel_lib angles CSV path (default: default.csv)')
    parser.add_argument('-s', '--sync', action="store_true", help='Generate kernels concurrently or not')

    parser.add_argument('-klib', '--kernel_folder', type=str, help='sgmap_multiple kernel folder path')
    parser.add_argument('-out', '--output_folder', type=str, help='sgmap_multiple output folder path')

    parser.add_argument('viirs_file', type=str, nargs='?', help='VIIRS image file path')

    parser.add_argument("--debug", help="increase output verbosity", action="store_true")

    args = parser.parse_args()
    if args.debug:
        logger.setLevel(logging.DEBUG)
    if args.action == 'gui':
        skyglow.main()
    else:
        if not args.viirs_file:
            raise ValueError('VIIRS file path is required for command line interface')

        if args.action == 'sgmap_single':
            zen, azi, krn_file = 0, 0, ''
            if not args.latitude:
                raise ValueError('Latitude is required for sgmap_single')
            if args.kernel_file:
                print('Using kernel file {} for sgmap_single'.format(args.kernel_file))
                krn_file = args.kernel_file
            darksky.sgmapper(args.latitude, args.clarity, zen, azi, args.viirs_file, krn_file)
        elif args.action == 'kernel_lib':
            if not args.latitude:
                raise ValueError('Latitude is required for kernel_lib')
            if not args.sync:
                static_args = [(args.latitude, args.clarity, args.viirs_file)]
                with open(args.angles_csv, "rb") as f:
                    angle_list = loadtxt(f, delimiter=",",skiprows=1)
                # cartesian product of arguments
                args_product = [(x[0], x[1], y[0], y[1], x[2]) for x in static_args for y in angle_list]
                p = Pool()
                try:
                    p.map_async(darksky.krn_unpack, args_product).get(9999999)
                except KeyboardInterrupt:
                    p.terminate()
            else:
                darksky.krnlibber(args.latitude, args.clarity, args.angles_csv, args.viirs_file)
        elif args.action == 'sgmap_multiple':
            if not args.kernel_folder:
                raise ValueError('Kernel folder path is required for sgmap_multiple')
            if not args.output_folder:
                raise ValueError('Output folder path is required for sgmap_multiple')
            darksky.multisgmapper(args.viirs_file, args.kernel_folder, args.output_folder)

if __name__ == '__main__':
    main()
