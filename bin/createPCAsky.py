#!/usr/bin/env python3

import argparse
import cubePCA
from astropy.io import fits as pyfits

def main():
    parser = argparse.ArgumentParser(description="Script to create a PCA library of spectra from an datacube given a mask")
    parser.add_argument('input_cube', metavar='CUBE_IN', type=str,
                        help='Input FITS datacube from which sky spectra are selected')
    parser.add_argument('PCA_out', metavar='PCA_OUT', type=str,
                        help='Output FITS file with PCA component spectra')
    parser.add_argument('sky_mask', metavar='MASK', type=str, help='sky region mask')
    parser.add_argument('-e', '--extension', type=int, default=1,
                        help='extension in which the data should be taken from')
    parser.add_argument('-s', '--spectra', type=int, default=20000,
                        help='Maximum random subset of spectra used for PCA analysis if number of selected spectra are larger than this number')
    parser.add_argument('-f', '--filter_width', type=int, default=50,
                        help='Size of median filter in wavelength direction to remove continuum signal before sky residual subtraction')
    parser.add_argument('--verbose', action='store_true', help='Set if infos are printed to the command line')

    args = parser.parse_args()


    if args.verbose:
        print('Opening data cube {} for processing '.format(args.input_cube))
    cube = cubePCA.IFUCube(args.input_cube,extension=args.extension)

    if args.verbose:
        print('Opening sky mask {} for processing '.format(args.sky_mask))
    mask = cubePCA.MASKimg(args.sky_mask)
    cube.getNANMask()

    if args.verbose:
        print('Creating PCA spectral libary. This may take a while... ')
    PCA_out,sky = cube.create_PCA_sky(mask, cont_filt=args.filter_width, spectra=args.spectra)
    hdu = pyfits.PrimaryHDU(PCA_out)
    hdu.writeto(args.PCA_out,overwrite=True)

main()