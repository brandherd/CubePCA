#!/usr/bin/env python3
import sys
import argparse
import cubePCA
from astropy.io import fits as pyfits
from tqdm import tqdm

__author__ = "Bernd Husemann"
__copyright__ = "Copyright 2020, Bernd Husemann"
__license__ = "MIT"
__version__ = "1.0"
__maintainer__ = "Bernd Husemann"
__email__ = "berndhusemannQgmx.de"
__status__ = "Production"


def main():
    parser = argparse.ArgumentParser(description="Script to subtract sky residuals from a datacube by creating a PCA spectral library")
    parser.add_argument('input_cube', metavar='CUBE_IN', type=str,  help='Input FITS datacube from which sky spectra are selected')
    parser.add_argument('output_cube', metavar='CUBE_OUT', type=str, help='Output FITS datacube with sky residuals subtracted')
    parser.add_argument('PCA_in', metavar='PCA_in', type=str, help='Output FITS file with PCA component spectra')
    parser.add_argument('-m', '--wave_mask', type=str, default='',help='File name of the wavelength selection for the PCA subtraction. If not given a default wavelength set is being used. An example is given in the wave.txt file.')
    parser.add_argument('-e', '--extension', type=int, default=1, help='extension in which the data should be taken from (Default: 1')
    parser.add_argument('-c', '--components', type=int, default=40, help='Number of PCA components to be used (Default: 40)')
    parser.add_argument('-f', '--filter_width', type=int, default=50, help='Size of median filter in wavelength direction to remove continuum signal before sky residual subtraction (Default: 50)')
    parser.add_argument('-p', '--processes', type=str, default='auto', help='number of processes to be used for multiprocessing, auto for maximum otherwise put number (Default: auto)')
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
        pbar = tqdm(total=cube.getSpax(), file=sys.stdout)
    else:
        pbar = None
    cube.subtract_PCA_sky(PCA, cont_filt=args.filter_width, components=args.components, file_wavemask= args.wave_mask, max_cpu=args.processes, pbar = pbar, verbose=args.verbose)

    if args.verbose:
        print('Store cleaned cube at {}.'.format(args.output_cube))
    cube.writeFits(args.output_cube)

    if args.verbose:
        print('Done')
if __name__ == '__main__':
    main()
