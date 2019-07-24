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

    def masked(self):
        dim = self.__mask[self.__extension].data.shape
        indices = numpy.indices(dim)
        indices_masked_y = indices[0][self.__mask[self.__extension].data.astype('bool')]
        indices_masked_x = indices[1][self.__mask[self.__extension].data.astype('bool')]
        return indices_masked_y,indices_masked_x



class IFUCube:
    def __init__(self, filename, extension=0):
        self.__hdu = pyfits.open(filename)
        self.extension = extension
        self.__header = self.__hdu[self.extension].header
        self.__dim = (self.__header['NAXIS3'], self.__header['NAXIS2'], self.__header['NAXIS1'])

    def subRSS(self,index_x,index_y):
        RSS = self.__hdu[self.extension].data[:,index_y,index_x]
        return RSS

    def writeFits(self,fileout):
        self.__hdu.writeto(fileout,overwrite=True)

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

        (indices_y,indices_x)= sky_mask.masked()
        indices = numpy.arange(len(indices_x))
        numpy.random.shuffle(indices)
        indices_y = indices_y[indices]
        indices_x = indices_x[indices]
        smax = numpy.min([spectra, len(indices_x)])
        sky = self.subRSS(indices_x[:smax],indices_y[:smax]).T
        sky_smooth = ndimage.filters.median_filter(sky,(1,cont_filt))
        sky = sky - sky_smooth
        M = numpy.dot(sky, sky.T).T
        e, EV = numpy.linalg.eigh(M)
        tmp = numpy.dot(sky.T, EV).T
        PCA_out = tmp[::-1]

        return PCA_out,sky

    def subtract_PCA_sky(self,PCA_sky,cont_filt=50,components=100,verbose=True):
        pca_specs = PCA_sky[:components, :] * 10000.0
        wave = self.getWave()
        select = ((wave>5565) & (wave<5585)) | ((wave>5875) & (wave<5905)) | ((wave>6290) & (wave<6310))| ((wave>6350) & (wave<6370))| ((wave>6450) & (wave<6610)) | ((wave>6800) & (wave<7020)) | ((wave>7220) & (wave<8160)) | ((wave>8250) & (wave<8680))| ((wave>8740) & (wave<9180))
        max_it = self.__dim[2]*self.__dim[1]

        bar_length = 30
        m=0
        for x in range(self.__dim[2]):
            for y in range(self.__dim[1]):
                spec = self.__hdu[self.extension].data[:,y,x]
                nan = numpy.isnan(spec)
                if numpy.sum(nan) < self.__dim[0]:
                    spec[nan] = 0
                    smooth_spec=ndimage.filters.median_filter(spec,(cont_filt))
                    out=numpy.linalg.lstsq(pca_specs[:,select].T,(spec-smooth_spec)[select])
                    spec_sky = numpy.dot(pca_specs[:, select].T, out[0])
                    self.__hdu[self.extension].data[select,y,x] = spec[select]-spec_sky
                m +=1
                if (m%100==0 or m== max_it-1) and verbose:
                    show_progress_bar(bar_length, m, max_it)
        print('\n')

    #def subtract_PCA_sky(self,PCA_sky):


#cube = IFUCube('/home/husemann/Projects/CARS/MUSE/HE0433-1028/HE0433-1028.unbinned.fits.gz')
#mask = MASKimg('/home/husemann/Projects/CARS/MUSE/HE0433-1028/SKY_MASK.fits')
#PCA_out,sky = cube.create_PCA_sky(mask)
#print("start sky subtraction")
#cube.subtract_PCA_sky(PCA_out)
#cube.writeFits('/home/husemann/Projects/CARS/MUSE/HE0433-1028/HE0433-1028.test.fits')

#def create_PCA_sky(cube_in,PCA_out,sky_mask,cont_filt=50,spectra=20000,parallel='auto'):

#    hdu = pyfits.open(cube_in)
#    cube = hdu[1].data
#    dim = cube.shape
#    hdu.close()

#    hdu = pyfits.open(sky_mask)
#    mask = hdu[0].data
    
#    nan = numpy.sum(numpy.isnan(cube[5:-5,:,:]),0)
#    cube[numpy.isnan(cube)] = 0.0
#    sky = cube[:,(mask==1) & (nan==0)].T
#    indices=numpy.arange(sky.shape[0])
#    numpy.random.shuffle(indices)
#    smax=numpy.min([spectra,sky.shape[0]])
#    out = (((sky)/10000.0).astype(numpy.float64))[indices[:smax],:]
#    if parallel=='auto':
#       	cpus = cpu_count()
#    else:
#      	cpus = int(parallel)
#    if cpus>1:
#    	pool = Pool(processes=cpus)
#    results=[]
#    for i in xrange(smax):
#	    results.append(pool.apply_async(ndimage.filters.median_filter,args=(out[i,:],cont_filt)))
#	pool.close()
#	pool.join()
#	smooth_out = numpy.zeros_like(out)
#	for c in xrange(smax):
#        spec=results[c].get()
#        smooth_out[c,:] = spec
#    else:
#	    smooth_out=ndimage.filters.median_filter(out,(1,cont_filt))
#    out = out-smooth_out
#    #hdu = pyfits.PrimaryHDU(out)
#    #hdu.writeto('test.fits',clobber=True)
#    M= numpy.dot(out,out.T).T
#    e,EV = numpy.linalg.eigh(M)
#    tmp = numpy.dot(out.T,EV).T
#    V=tmp[::-1]
    #S=numpy.sqrt(e)[::-1]

#    hdu = pyfits.PrimaryHDU(V[:dim[0],:])
#    hdu.writeto(PCA_out,clobber=True)

#def fit_PCA(cube,idx_x,idx_y,pca_specs,cont_filt,select):
#    dim = cube.shape
#    clean_cube = numpy.zeros((numpy.sum(select),dim[1]))
#    for m in xrange(len(idx_x)):
        #print m
#        spec = cube[:,m]
#        if numpy.sum(numpy.isnan(spec))<10:
#            spec[numpy.isnan(spec)] = 0
#            smooth_spec=ndimage.filters.median_filter(spec,(cont_filt))
#            out=numpy.linalg.lstsq(pca_specs[:,select].T,(spec-smooth_spec)[select])
#            spec_sky = numpy.dot(pca_specs[:,select].T,out[0])
#        else:
#            spec_sky = 0.0
#        clean_cube[:,m] = spec[select]-spec_sky
        #if m>20000:
        #    break
#    return clean_cube

#def subtract_PCA_sky(cube_in,cube_out,PCA_spec,components=150,cont_filt=40,parallel='auto'):
#    hdu = pyfits.open(PCA_spec)
#    pca_specs = hdu[0].data[:components,:]*10000.0
#    hdu.close()
#    hdu = pyfits.open(cube_in)
#    cube = hdu[1].data
#    hdr = hdu[1].header
#    dim = cube.shape
#    wave = numpy.arange(dim[0])*hdr['CD3_3']+hdr['CRVAL3']
#    select = ((wave>5565) & (wave<5585)) | ((wave>5875) & (wave<5905)) | ((wave>6290) & (wave<6310))| ((wave>6350) & (wave<6370))| ((wave>6450) & (wave<6610)) | ((wave>6800) & (wave<7020)) | ((wave>7220) & (wave<8160)) | ((wave>8250) & (wave<8680))| ((wave>8740) & (wave<9180))

#    clean_cube = numpy.zeros((numpy.sum(select),dim[1],dim[2]))
#    clean_cube[:,:,:] = cube[select,:,:]
#    if parallel=='auto':
#	    cpus = cpu_count()
#    else:
#            cpus = int(parallel)
#    idx = numpy.indices((dim[1],dim[2]))
#    idx_x = idx[1].flatten()
#    idx_y = idx[0].flatten()
#    if cpus>1:
#	    pool = Pool(processes=cpus)
#	    idx_y_split = numpy.array_split(idx_y,cpus)
#            idx_x_split = numpy.array_split(idx_x,cpus)
#            results=[]
#            for c in xrange(cpus):
#                results.append(pool.apply_async(fit_PCA,args=(cube[:,idx_y_split[c],idx_x_split[c]],idx_x_split[c],idx_y_split[c],pca_specs,cont_filt,select)))
#	    pool.close()
#            pool.join()
#            for c in xrange(cpus):
#                out=results[c].get()
#                clean_cube[:,idx_y_split[c],idx_x_split[c]] = out[:,:]
#    else:

#       clean_rss = fit_PCA(cube[:,idx_y,idx_x],idx_y,idx_x,pca_specs,cont_filt,select)
#       clean_cube[:,idx_y,idx_x] = clean_rss[:,:]

#    cube[select,:,:] = clean_cube
#    hdu.writeto(cube_out,clobber=True)



#def _pickle_method(method):
#      func_name = method.im_func.__name__
#      obj = method.im_self
#      cls = method.im_class
#      if func_name.startswith('__') and not func_name.endswith('__'):
#	cls_name = cls.__name__.lstrip('_')
#	if cls_name: func_name = '_' + cls_name + func_name
#      return _unpickle_method, (func_name, obj, cls)
    
#def _unpickle_method(func_name, obj, cls):
#      for cls in cls.mro():
#	try:
#	  func = cls.__dict__[func_name]
#	except KeyError:
#	  pass
#	else:
#	  break
#     return func.__get__(obj, cls)
#copy_reg.pickle(MethodType,_pickle_method, _unpickle_method)
