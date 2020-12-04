from astropy.io import fits as pyfits
import math
import numpy
from scipy import ndimage


def show_progress_bar(bar_length, completed, total):
    bar_length_unit_value = (total / bar_length)
    completed_bar_part = math.ceil(completed / bar_length_unit_value)
    progress = "*" * completed_bar_part
    remaining = " " * (bar_length - completed_bar_part)
    percent_done = "%.2f" % ((completed / total) * 100)
    print(f'[{progress}{remaining}] {percent_done}%', end='\r')


class MASKimg:
    def __init__(self, filename):
        self.__mask = pyfits.open(filename)
        if self.__mask[0].header['NAXIS']==0:
            self.__extension = 1
        else:
            self.__extension = 0

    def masked(self,nan_mask = None):
        dim = self.__mask[self.__extension].data.shape
        indices = numpy.indices(dim)
        if nan_mask is None:
            indices_masked_y = indices[0][self.__mask[self.__extension].data.astype('bool')]
            indices_masked_x = indices[1][self.__mask[self.__extension].data.astype('bool')]
        else:
            indices_masked_y = indices[0][(self.__mask[self.__extension].data.astype('bool')) & (nan_mask==False)]
            indices_masked_x = indices[1][(self.__mask[self.__extension].data.astype('bool')) & (nan_mask==False)]
        return indices_masked_y,indices_masked_x

class WAVEmask:
    def __init__(self, filename=''):
        if filename!='':
            file_in = open(filename,'r')
            lines = file_in.readlines()
            wave_start = []
            wave_end = []
            for i in range(len(lines)):
                line = lines[i].split()
                if len(line)==2:
                    wave_start.append(float(line[0]))
                    wave_end.append(float(line[1]))
                else:
                    break
        else:
            wave_start = [5565 ,5875, 6290, 6350, 6450, 6800, 7220, 8250, 8740]
            wave_end =[5585, 5905, 6310, 6370, 6610, 7020, 8160, 8680, 9180]

        self.__wave_start = numpy.array(wave_start)
        self.__wave_end = numpy.array(wave_end)

    def mask(self, wave):
        wave_mask = numpy.zeros(len(wave),dtype='bool')
        for i in range(len(self.__wave_start)):
            wave_mask = wave_mask | ((wave>=self.__wave_start[i]) & (wave<=self.__wave_end[i]))
        return wave_mask


class IFUCube:
    def __init__(self, filename, extension=0):
        self.__hdu = pyfits.open(filename)
        self.extension = extension
        self.__header = self.__hdu[self.extension].header
        self.__dim = (self.__header['NAXIS3'], self.__header['NAXIS2'], self.__header['NAXIS1'])
        self.__badmask = None

    def subRSS(self,index_x,index_y):
        RSS = self.__hdu[self.extension].data[:,index_y,index_x]
        return RSS

    def writeFits(self,fileout):
        self.__hdu.writeto(fileout,overwrite=True)

    def getNANMask(self,min_slice=5, max_slice=20):
        self.__badmask = numpy.sum(numpy.isnan(self.__hdu[self.extension].data[5:20,:,:]),0)>0.0

    def replaceNAN(self,value=0.0):
        select_nan = numpy.isnan(self.__hdu[self.extension].data)
        self.__hdu[self.extension].data[select_nan] = value

    def getWave(self):
        try:
            crpix = self.__header['CRPIX3']
        except KeyError:
            crpix = 1
        crval = self.__header['CRVAL3']
        try:
            cdelt = self.__header['CD3_3']
        except KeyError:
            cdelt = self.__header['CDELT3']
        wave = (numpy.arange(self.__dim[0])-(crpix-1))*cdelt+crval
        return wave


    def create_PCA_sky(self, sky_mask, cont_filt=50, spectra=20000, parallel='auto'):

        (indices_y,indices_x)= sky_mask.masked(nan_mask = self.__badmask)
        indices = numpy.arange(len(indices_x))
        numpy.random.shuffle(indices)
        indices_y = indices_y[indices]
        indices_x = indices_x[indices]
        smax = numpy.min([spectra, len(indices_x)])
        self.replaceNAN()
        sky = self.subRSS(indices_x[:smax],indices_y[:smax]).T
        sky_smooth = ndimage.filters.median_filter(sky,(1,cont_filt))
        sky = sky - sky_smooth
        M = numpy.dot(sky, sky.T).T
        e, EV = numpy.linalg.eigh(M)
        tmp = numpy.dot(sky.T, EV).T
        PCA_out = tmp[::-1]

        return PCA_out,sky

    def subtract_PCA_sky(self, PCA_sky, cont_filt=50, components=100, file_wavemask='', verbose=True):
        pca_specs = PCA_sky[:components, :] * 10000.0
        wave_mask = WAVEmask(file_wavemask)

        select_wave = wave_mask.mask(self.getWave())
        max_it = self.__dim[2]*self.__dim[1]

        bar_length = 30
        m=0
        for x in range(self.__dim[2]):
            for y in range(self.__dim[1]):
                spec = self.__hdu[self.extension].data[:,y,x]
                nan = numpy.isnan(spec)
                spec[nan] = 0
                select = select_wave
                print(numpy.sum(select))
                print((pca_specs[:,select].T))
                print (spec[select])
                print(cont_filt)
                smooth_spec=ndimage.filters.median_filter(spec,(cont_filt))
                print(smooth_spec)
                out=numpy.linalg.lstsq(pca_specs[:,select].T,(spec-smooth_spec)[select])
                spec_sky = numpy.dot(pca_specs[:, select].T, out[0])
                self.__hdu[self.extension].data[select,y,x] = spec[select]-spec_sky
                m +=1
                if (m%100==0 or m== max_it-1) and verbose:
                    show_progress_bar(bar_length, m, max_it)
        print('\n')