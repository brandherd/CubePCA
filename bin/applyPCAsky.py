#!/usr/bin/env python3

import argparse
from cubePCA import cubePCA
from astropy.io import fits as pyfits


def main():
    parser = argparse.ArgumentParser(description="Script to subtract sky residuals from a datacube by creating a PCA spectral library")
    parser.add_argument('input_cube', metavar='CUBE_IN', type=str,  help='Input FITS datacube from which sky spectra are selected')
    parser.add_argument('output_cube', metavar='CUBE_OUT', type=str, help='Output FITS datacube with sky residuals subtracted')
    parser.add_argument('PCA_in', metavar='PCA_in', type=str, help='Output FITS file with PCA component spectra')
    parser.add_argument('-m', '--wave_mask', type=str, default='',help='File name of the wavelength selection for the PCA subtraction. If not given a default wavelength set is being used. An example is given in the wave.txt file.')
    parser.add_argument('-e', '--extension', type=int, default=1, help='extension in which the data should be taken from')
    parser.add_argument('-c', '--components', type=int, default=100, help='Number of PCA components to be used')
    parser.add_argument('-f', '--filter_width', type=int, default=50, help='Size of median filter in wavelength direction to remove continuum signal before sky residual subtraction')
    parser.add_argument('--verbose', action='store_true', help='Set if infos are printed to the command line')

    args = parser.parse_args()
    if args.verbose:
        print('Opening data cube {} for processing '.format(args.input_cube))
    cube = cubePCA.IFUCube(args.input_cube,extension=args.extension)


    if args.verbose:
        print('Reading PCA spectral libary from disc')
    PCA = pyfits.open(args.PCA_in)[0].data

    if args.verbose:
        print('Start subtracting sky line residuals.')
    cube.subtract_PCA_sky(PCA, cont_filt=args.filter_width, components=args.components, file_wavemask= args.wave_mask, verbose=args.verbose)

    if args.verbose:
        print('Store cleaned cube at {}.'.format(args.output_cube))
    cube.writeFits(args.output_cube)

    if args.verbose:
        print('Done')

main()
