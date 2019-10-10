from astropy.io import fits


def hdr3Dto2D(hdr3D,verbose=True):
	"""
	Removing the wavelength component of a header, i.e., converting
	the hdr from 3D (lambda, dec, ra) to 2D (dec, ra)

	--- INPUT ---
	hdr3d       The header object to convert from (lambda, dec, ra) to (dec, ra)
	verbose     Toggle verbosity

	"""
	for key in hdr3D.keys():
		if any(a in key for a in ('S3', '_3')):
			del hdr3D[key]
			
	return hdr3D


filename = 'dh_SL1_cube'
wvlnum = 5 # 0-indexed
with fits.open(filename+'.fits') as hdul0:
	data0 = hdul0[0].data[wvlnum]
	hdr0 = hdul0[0].header
	data1 = hdul0[1].data[0][0] # 0D list; data type is FITS_rec([([[*],[*],...,])], dtype=...)
	## WCS keywords
#	CTYPE1 = hdr0['CTYPE1']
#	CTYPE2 = hdr0['CTYPE2']
#	CRPIX1 = hdr0['CRPIX1'] # ref ax
#	CRPIX2 = hdr0['CRPIX2'] # ref ay
#	CRVAL1 = hdr0['CRVAL1'] # ref ra
#	CRVAL2 = hdr0['CRVAL2'] # ref dec
#	CDELT1 = hdr0['CDELT1'] # delta ra
#	CDELT2 = hdr0['CDELT2'] # delta dec
#	PC1_1 = hdr0['PC1_1'] # rot matrix
#	PC2_1 = hdr0['PC2_1']
#	PC1_2 = hdr0['PC1_2']
#	PC2_2 = hdr0['PC2_2']
	## alternative header extracting
	hdr = hdr3Dto2D(hdr0)
	wvl = data1[0] + (data1[1] - data1[0]) * wvlnum
	wvl = wvl[0]
	hdr['CRPIX3'] = (wvlnum, "sliced wavelength number")
	hdr['CRVAL3'] = (wvl, "current wavelength is {} um".format(wvl))
	hdr.add_comment("The cube is sliced at {} um.".format(wvl))
## new fits file creating
primary_hdu = fits.PrimaryHDU(data0)
hdul = fits.HDUList(primary_hdu)
hdul[0].header = hdr # copy all (2D) header info
#hdr = hdul[0].header
#hdr['CTYPE1'] = CTYPE1
#hdr['CTYPE2'] = CTYPE2
#hdr['CRPIX1'] = CRPIX1
#hdr['CRPIX2'] = CRPIX2
#hdr['CRVAL1'] = CRVAL1
#hdr['CRVAL2'] = CRVAL2
#hdr['CDELT1'] = CDELT1
#hdr['CDELT2'] = CDELT2
#hdr['PC1_1'] = PC1_1
#hdr['PC2_1'] = PC2_1
#hdr['PC1_2'] = PC1_2
#hdr['PC2_2'] = PC2_2
#hdr['EQUINOX'] = 2000.0
hdul.writeto(filename+'_s.fits', output_verify='ignore', overwrite=True)
