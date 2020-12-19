#!/usr/bin/env python3
import os
import numpy
import argparse
from astropy.io import fits as pyfits
from matplotlib import pyplot as plt
import matplotlib

__author__ = "Bernd Husemann"
__copyright__ = "Copyright 2020, Bernd Husemann"
__license__ = "MIT"
__version__ = "1.0"
__maintainer__ = "Bernd Husemann"
__email__ = "berndhusemann@gmx.de"
__status__ = "Production"

class maskArray:
    def __init__(self, dim):
        self.mask = numpy.zeros(dim,dtype=numpy.uint16)
        self.dim = dim
        self.temp_corners_x = []
        self.temp_corners_y = []

    def replaceMask(self, mask):
        self.mask = mask
        self.dim = mask.shape

    def addRegion(self,xpos,ypos):
        self.mask[min(ypos):max(ypos)+1,min(xpos):max(xpos)+1] = 1

    def removeRegion(self,xpos,ypos):
        self.mask[min(ypos):max(ypos)+1,min(xpos):max(xpos)+1] = 0

    def loadMask(self,filename):
        maskhdu = pyfits.open(filename)
        data = maskhdu[0].data
        if self.dim==data.shape:
            self.mask = data
        else:
            raise ImportError('Existing mask is not matching dimensions of reference data cube')


    def saveMask(self,filename):
        maskHDU = pyfits.PrimaryHDU(self.mask)
        maskHDU.writeto(filename,overwrite=True)

    def plotMask(self,image):
        self.fig = plt.figure()
        self.vmin=0
        self.vmax=numpy.nanmax(image)/10.0
        self.xmin=0
        self.xmax=self.dim[1]
        self.ymin = 0
        self.ymax = self.dim[0]
        self.ax = self.fig.add_subplot(111)
        self.ax.set_title("press left -> mask; press right -> unmask; press c -> change cuts")
        self.im = self.ax.imshow(image, origin='lower', cmap=matplotlib.cm.viridis, interpolation='nearest', vmin=self.vmin, vmax=self.vmax, zorder=1)
        self.contour = self.ax.contourf(self.mask,levels=[0.9,1.1],colors='r',alpha=0.5,zorder=2)
        self.ax.set_xlabel('x pixel')
        self.ax.set_ylabel('y pixel')
        self.fig.tight_layout()
        self._connect()
        plt.show()

    def updatePlot(self):
        for tp in self.contour.collections:
            tp.remove()
        self.contour = self.ax.contourf(self.mask, levels=[0.9, 1.1], colors='r', alpha=0.5, zorder=2)
        self.im.set_clim(vmin=self.vmin,vmax=self.vmax)
        self.ax.set_xlim([self.xmin,self.xmax])
        self.ax.set_ylim([self.xmin, self.ymax])
        plt.draw()

    def press_on(self,event):
        if event.inaxes != self.ax.axes:
            print('Outside of image region')
            return
        else:
            self.temp_corners_x.append(int(round(event.xdata)))
            self.temp_corners_y.append(int(round(event.ydata)))

    def press_off(self,event):
        if event.inaxes != self.ax.axes:
            print('Outside of image region')
            self.temp_corners_x = []
            self.temp_corners_y = []
            return
        else:
            self.temp_corners_x.append(int(round(event.xdata)))
            self.temp_corners_y.append(int(round(event.ydata)))
        if len(self.temp_corners_x) == 2:
            if event.button == 1:
                self.addRegion(self.temp_corners_x, self.temp_corners_y)
            elif event.button == 3:
                self.removeRegion(self.temp_corners_x, self.temp_corners_y)
            self.updatePlot()
            self.temp_corners_x = []
            self.temp_corners_y = []
        else:
            self.temp_corners_x = []
            self.temp_corners_y = []

    def on_key(self,event):
        if event.key == 'c':
            print('Current cut levels are: %f,%f' %(self.vmin, self.vmax))
            new_limits = input('Enter new cut levels: ')
            self.vmin = float(new_limits.split(',')[0])
            self.vmax = float(new_limits.split(',')[1])
            self.updatePlot()
        if event.key == 'x':
            print('Current x-axis visibility is: %d,%d' %(self.xmin, self.xmax))
            new_limits = input('Enter new x visiblity levels or enter to unzoom: ')
            if len(new_limits.split(','))==2:
                self.xmin = int(new_limits.split(',')[0])
                self.xmax = int(new_limits.split(',')[1])
            else:
                self.xmin = 0
                self.xmax = self.dim[1]
            self.updatePlot()
        if event.key == 'y':
            print('Current y-axis visibility is: %d,%d' %(self.ymin, self.ymax))
            new_limits = input('Enter new y visiblity levels or enter to unzoom: ')
            if len(new_limits.split(','))==2:
                self.ymin = int(new_limits.split(',')[0])
                self.ymax = int(new_limits.split(',')[1])
            else:
                self.ymin = 0
                self.ymax = self.dim[0]
            self.updatePlot()

    def _connect(self):
        self.ax.figure.canvas.mpl_connect('button_press_event', self.press_on)
        self.ax.figure.canvas.mpl_connect('button_release_event', self.press_off)
        self.fig.canvas.mpl_connect('key_press_event', self.on_key)



def main():
    parser = argparse.ArgumentParser(description="Script to visually create a mask file from a data cube. " 
                                                 "Keep left or right mouse button pressed while drawing a rectangle from one to the other corner. "
                                                 "The left button will add a region to the mask and the right button will remove a region from the mask."
                                                 "Press c button to change the linear intensity scaling of the image."
                                                 "Press x button to set a different x-axis visiblity for plotting."
                                                 "Press y button to set a different y-axis visiblity for plotting.")
    parser.add_argument('input_cube', metavar='CUBE_IN', type=str,
                        help='Input FITS datacube for which a mask file will be generated')
    parser.add_argument('mask_out', metavar='MASK_OUT', type=str,
                        help='Output FITS file with mask')
    parser.add_argument('-e', '--extension', type=int, default=1,
                        help='extension in which the data should be taken from  data cube (Default: 1)')
    parser.add_argument('-sw', '--startw', type=int, default=0,
                        help='Start wavelength pixel used for collapsing the cube to an image. (Default: 0)')
    parser.add_argument('-ew', '--endw', type=int, default=-1,
                        help='End wavelength pixel used for collapsing the cube to an image. (Default: -1)')
    parser.add_argument('-m', '--method', type=str, default='median',
                        help='Method to collapse the data along the wavelength axis within the given axis range. Available: median, mean, sum, max. (Default: median)')



    args = parser.parse_args()
    extension = args.extension
    method = args.method
    hdu = pyfits.open(args.input_cube)
    if extension is not None:
        data = hdu[extension].data
    else:
        if len(hdu[0].data)>1:
            data = hdu[0].data
        else:
            data = hdu[1].data
    if len(data)<2:
        raise ImportError('Data not in suitable format, check fits file or provide correct data extension')
    wmin=args.startw
    wmax=args.endw
    if method=='median':
        collapseed_image = numpy.nanmedian(data[wmin:wmax,:,:],axis=0)
    elif method=='sum':
        collapseed_image = numpy.nansum(data[wmin:wmax, :, :], axis=0)
    elif method=='mean':
        collapseed_image = numpy.nanmean(data[wmin:wmax, :, :], axis=0)
    elif method=='max':
        collapseed_image = numpy.nanmax(data[wmin:wmax, :, :], axis=0)
    else:
        raise RuntimeError('Collapsing method %s is not defined. Please use either median, mean, sum or max.'%(method))
    dim_data = data.shape
    mask = maskArray(dim_data[1:])
    if os.path.isfile(args.mask_out):
        mask.loadMask(args.mask_out)
    mask.plotMask(collapseed_image)
    mask.saveMask(args.mask_out)


if __name__ == '__main__':
    main()











