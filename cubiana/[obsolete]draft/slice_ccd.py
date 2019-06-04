from astropy.nddata import CCDData

wvlnum = 5
filename = 'dh_SL2_cube'
data = CCDData.read(filename+'.fits')
print("---------- Ignore the info above ----------")
ccd = CCDData(data, mask=None, unit=None)
ccd1 = ccd[wvlnum,:,:]
print("The {}th wavelength map is extracted.".format(wvlnum))
ccd1.write(filename+'_s.fits', overwrite=True)
