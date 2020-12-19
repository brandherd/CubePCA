createPCAsky.py MUSE_example.fits.gz PCA_sky.fits SKY_mask.fits -e 0 --verbose
applyPCAsky.py MUSE_example.fits.gz MUSE_example_clean.fits PCA_sky.fits -e 0 -c 20 -m wave_mask.txt -p 4 --verbose
