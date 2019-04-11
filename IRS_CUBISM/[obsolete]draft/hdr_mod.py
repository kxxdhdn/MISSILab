from astropy.io import fits


filename = 'sh_SL2_cube'
with fits.open(filename+'.fits') as hdul0:
	data0 = hdul0[0].data
	hdr0 = hdul0[0].header
	NAXIS3 = hdr0['NAXIS3']
	hdr0['EQUINOX'] = 2000.0
	table = hdul0[1].data
	wvl = []
	for i in range(NAXIS3):
		wvl.append(table[0][0][i][0])
	hdr1 = hdul0[1].header
primary_hdu = fits.PrimaryHDU(header=hdr0, data=data0)
hdul = fits.HDUList(primary_hdu)
hdu = fits.ImageHDU(data=wvl, header=hdr1, name='wvl') # add table
hdul.append(hdu)
print(hdul[0].header['EQUINOX'])
hdul.writeto(filename+'_m.fits', output_verify='ignore', overwrite=True)
