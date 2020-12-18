# CubePCAsky Package

This is a simple package for cleaning astronomical IFU osbservations from sky line residuals or any systematic artefacts. 
It has been developed mainly for VLT/MUSE optical IFU data, but may also be applicable to other data sets. 
The software provides 3 command-line tools to performs the necessary tasks:
* `createPCAsky.py` to create PCA component spectra
* `applyPCAsky.py` to clean a cube using the PCA component spectra
* `subtractPCAsky.py` which combines both step into a single command

A common sequence of commands would be as follows 
1. `createPCAsky.py DATACUBE_IN.fits PCA_LIBRARY_OUT.fits SKY_MASK.fits --verbose` 
2. `applyPCAsky.py DATACUBE_IN.fits DATACUBE_CLEAN.fits PCA_LIBRARY_OUT.fits -m WAVE_MASK.txt --verbose`

There are a few more parameters available for each command to control the process. Please use `-h` option for more details. 
In most cases the default parameter should be sufficient at least in case of VLT/MUSE observations. 

A SKY_MASK.fits file has to be created either manually or automatically which needs to have the same
spatial dimensions as the input data cube. Regions to be considered as pure sky for the PCA analysis 
need to set to 1 and others to 0. 




