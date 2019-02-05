from astropy.io import fits
from astropy.wcs import WCS


def rhdr(fits_path, fits_filename):
	with fits.open(fits_path+fits_filename+'.fits') as hdul:
		##header
		hdr = hdul[0].header
		CTYPE1 = hdr['CTYPE1']
		CTYPE2 = hdr['CTYPE2']
		NAXIS1 = hdr['NAXIS1'] #nb ax
		NAXIS2 = hdr['NAXIS2'] #nb ay
		CRPIX1 = hdr['CRPIX1'] #ref ax
		CRPIX2 = hdr['CRPIX2'] #ref ay
		CRVAL1 = hdr['CRVAL1'] #ref ra
		CRVAL2 = hdr['CRVAL2'] #ref dec
		CDELT1 = hdr['CDELT1'] #delta ra
		CDELT2 = hdr['CDELT2'] #delta dec
		PC1_1 = hdr['PC1_1'] #rot matrix
		PC2_1 = hdr['PC2_1']
		PC1_2 = hdr['PC1_2']
		PC2_2 = hdr['PC2_2']

	w = WCS(hdr)
	return hdr, CTYPE1, CTYPE2, NAXIS1, NAXIS2, \
		CRPIX1, CRPIX2, CRVAL1, CRPIX2, CDELT1, CDELT2, \
		PC1_1, PC2_1, PC1_2, PC2_2, w
	##(5, 6, 5) total = 16



def rhdr3d(fits_path, fits_filename):
	with fits.open(fits_path+fits_filename+'.fits') as hdul:
		##header
		hdr = hdul[0].header
		CTYPE1 = hdr['CTYPE1']
		CTYPE2 = hdr['CTYPE2']
		NAXIS1 = hdr['NAXIS1'] #nb ax
		NAXIS2 = hdr['NAXIS2'] #nb ay
		NAXIS3 = hdr['NAXIS3'] #nb az/ wvl
		CRPIX1 = hdr['CRPIX1'] #ref ax
		CRPIX2 = hdr['CRPIX2'] #ref ay
		CRVAL1 = hdr['CRVAL1'] #ref ra
		CRVAL2 = hdr['CRVAL2'] #ref dec
		CDELT1 = hdr['CDELT1'] #delta ra
		CDELT2 = hdr['CDELT2'] #delta dec
		#PC1_1 = hdr['PC1_1'] #rot matrix
		#PC2_1 = hdr['PC2_1']
		#PC1_2 = hdr['PC1_2']
		#PC2_2 = hdr['PC2_2']

	#w = WCS(hdr).celestial
	return hdr, CTYPE1, CTYPE2, NAXIS1, NAXIS2, NAXIS3, \
		CRPIX1, CRPIX2, CRVAL1, CRPIX2, CDELT1, CDELT2, \
		#PC1_1, PC2_1, PC1_2, PC2_2, w
	##(6, 6, 5) total = 17



def cphdr(data, fits_path, fits_filename, \
	new_fits_path, new_fits_filename):
	with fits.open(fits_path+fits_filename+'.fits') as hdul:
		##load header
		hdr = hdul[0].header
		CTYPE1 = hdr['CTYPE1']
		CTYPE2 = hdr['CTYPE2']
		CRPIX1 = hdr['CRPIX1'] #ref ax
		CRPIX2 = hdr['CRPIX2'] #ref ay
		CRVAL1 = hdr['CRVAL1'] #ref ra
		CRVAL2 = hdr['CRVAL2'] #ref dec
		CDELT1 = hdr['CDELT1'] #delta ra
		CDELT2 = hdr['CDELT2'] #delta dec
		PC1_1 = hdr['PC1_1'] #rot matrix
		PC2_1 = hdr['PC2_1']
		PC1_2 = hdr['PC1_2']
		PC2_2 = hdr['PC2_2']
	##create new hdr
	primary_hdu = fits.PrimaryHDU(data)
	hdul = fits.HDUList(primary_hdu)
	hdr = hdul[0].header
	hdr['CTYPE1'] = CTYPE1
	hdr['CTYPE2'] = CTYPE2
	hdr['CRPIX1'] = CRPIX1
	hdr['CRPIX2'] = CRPIX2
	hdr['CRVAL1'] = CRVAL1
	hdr['CRVAL2'] = CRVAL2
	hdr['CDELT1'] = CDELT1
	hdr['CDELT2'] = CDELT2
	hdr['PC1_1'] = PC1_1
	hdr['PC2_1'] = PC2_1
	hdr['PC1_2'] = PC1_2
	hdr['PC2_2'] = PC2_2
	hdr['EQUINOX'] = 2000.0
	##write new fits
	hdul.writeto(new_fits_path+new_fits_filename+'.fits', overwrite=True)



def cphdr3d(data3d, table, tablename, comment, fits_path, fits_filename, \
	new_fits_path, new_fits_filename):
	with fits.open(fits_path+fits_filename+'.fits') as hdul:
		##load header
		hdr = hdul[0].header
		CTYPE1 = hdr['CTYPE1']
		CTYPE2 = hdr['CTYPE2']
		CRPIX1 = hdr['CRPIX1'] #ref ax
		CRPIX2 = hdr['CRPIX2'] #ref ay
		CRVAL1 = hdr['CRVAL1'] #ref ra
		CRVAL2 = hdr['CRVAL2'] #ref dec
		CDELT1 = hdr['CDELT1'] #delta ra
		CDELT2 = hdr['CDELT2'] #delta dec
		PC1_1 = hdr['PC1_1'] #rot matrix
		PC2_1 = hdr['PC2_1']
		PC1_2 = hdr['PC1_2']
		PC2_2 = hdr['PC2_2']
	##create new hdr
	hdr = fits.Header()
	hdr['OBSERVER'] = 'Spitzer IRS'
	hdr['COMMENT'] = comment
	hdr['CTYPE1'] = CTYPE1
	hdr['CTYPE2'] = CTYPE2
	hdr['CRPIX1'] = CRPIX1
	hdr['CRPIX2'] = CRPIX2
	hdr['CRVAL1'] = CRVAL1
	hdr['CRVAL2'] = CRVAL2
	hdr['CDELT1'] = CDELT1
	hdr['CDELT2'] = CDELT2
	hdr['PC1_1'] = PC1_1
	hdr['PC2_1'] = PC2_1
	hdr['PC1_2'] = PC1_2
	hdr['PC2_2'] = PC2_2
	hdr['EQUINOX'] = 2000.0
	primary_hdu = fits.PrimaryHDU(header=hdr, data=data3d)
	hdul = fits.HDUList(primary_hdu)
	##add table
	hdu = fits.ImageHDU(data=table, name=tablename)
	hdul.append(hdu)
	##write new fits
	hdul.writeto(new_fits_path+new_fits_filename+'.fits', overwrite=True)



def cd2pc(data3d, fits_path, fits_filename, \
	new_fits_path, new_fits_filename):
	with fits.open(fits_path+fits_filename+'.fits') as hdul:
		hdr = hdul[0].header
		CRPIX1 = hdr['CRPIX1']
		CRPIX2 = hdr['CRPIX2']
		CDELT1 = 1.45/3600.#hdr['CDELT1'] faut in the hdr
		CDELT2 = 1.45/3600.#hdr['CDELT2'] faut in the hdr
		CUNIT1 = hdr['CUNIT1']
		CUNIT2 = hdr['CUNIT2']
		CTYPE1 = hdr['CTYPE1']
		CTYPE2 = hdr['CTYPE2']
		CRVAL1 = hdr['CRVAL1']
		CRVAL2 = hdr['CRVAL2']
		CRPIX3 = hdr['CRPIX3']
		CRVAL3 = hdr['CRVAL3']
		CDELT3 = hdr['CDELT3']
		CTYPE3 = hdr['CTYPE3']
		BUNIT = hdr['BUNIT']
		OBJECT= hdr['OBJECT']
		LONPOLE = hdr['LONPOLE']
		LATPOLE = hdr['LATPOLE']
		RADESYS = hdr['RADESYS']
		EQUINOX = hdr['EQUINOX']
		CD1_1 = hdr['CD1_1']
		CD2_1 = hdr['CD2_1']
		CD1_2 = hdr['CD1_2']
		CD2_2 = hdr['CD2_2']
		CD3_3 = hdr['CD3_3']
		EXTNAME = hdr['EXTNAME']

	primary_hdu = fits.PrimaryHDU(data3d)
	hdul = fits.HDUList(primary_hdu)
	hdr = hdul[0].header
	hdr['CRPIX1'] = CRPIX1
	hdr['CRPIX2'] = CRPIX2
	hdr['CDELT1'] = CDELT1
	hdr['CDELT2'] = CDELT2
	hdr['CUNIT1'] = CUNIT1
	hdr['CUNIT2'] = CUNIT2
	hdr['CTYPE1'] = CTYPE1
	hdr['CTYPE2'] = CTYPE2
	hdr['CRVAL1'] = CRVAL1
	hdr['CRVAL2'] = CRVAL2
	#hdr['CRPIX3'] = CRPIX3
	#hdr['CRVAL3'] = CRVAL3
	#hdr['CDELT3'] = CDELT3
	#hdr['CTYPE3'] = CTYPE3
	hdr['BUNIT'] = BUNIT
	hdr['OBJECT'] = OBJECT
	hdr['LONPOLE'] = LONPOLE
	hdr['LATPOLE'] = LATPOLE
	hdr['RADESYS'] = RADESYS
	hdr['EQUINOX'] = EQUINOX
	hdr['PC1_1'] = CD1_1 / CDELT1
	hdr['PC2_1'] = CD2_1 / CDELT1
	hdr['PC1_2'] = CD1_2 / CDELT2
	hdr['PC2_2'] = CD2_2 / CDELT2
	hdr['EXTNAME'] = EXTNAME
	##write new fits
	hdul.writeto(new_fits_path+new_fits_filename+'.fits', overwrite=True)



def zoom(data, fits_path, fits_filename, \
	new_fits_path, new_fits_filename, \
	xmin, xmax, ymin, ymax, w):
	with fits.open(fits_path+fits_filename+'.fits') as hdul:
		##header
		hdr = hdul[0].header
		CTYPE1 = hdr['CTYPE1']
		CTYPE2 = hdr['CTYPE2']
		CRPIX1 = hdr['CRPIX1'] #ref ax
		CRPIX2 = hdr['CRPIX2'] #ref ay
		CRVAL1 = hdr['CRVAL1'] #ref ra
		CRVAL2 = hdr['CRVAL2'] #ref dec
		CDELT1 = hdr['CDELT1'] #delta ra
		CDELT2 = hdr['CDELT2'] #delta dec
		PC1_1 = hdr['PC1_1'] #rot matrix
		PC2_1 = hdr['PC2_1']
		PC1_2 = hdr['PC1_2']
		PC2_2 = hdr['PC2_2']
	##modification of coord ref
	nrpix1 = (xmax - xmin) / 2 #new CRPIXi (center by default)
	nrpix2 = (ymax - ymin) / 2
	ax = (xmax + xmin) / 2 #new CRPIXi in old coord
	ay = (ymax + ymin) / 2
	nrval1, nrval2 = w.all_pix2world(ax, ay, 1) #new CRVALi (w related to old hdr)
	##create new hdr
	primary_hdu = fits.PrimaryHDU(data)
	hdul = fits.HDUList(primary_hdu)
	hdr = hdul[0].header
	hdr['CTYPE1'] = CTYPE1
	hdr['CTYPE2'] = CTYPE2
	hdr['CRPIX1'] = nrpix1
	hdr['CRPIX2'] = nrpix2
	hdr['CRVAL1'] = float(nrval1)
	hdr['CRVAL2'] = float(nrval2)
	hdr['CDELT1'] = CDELT1
	hdr['CDELT2'] = CDELT2
	hdr['PC1_1'] = PC1_1
	hdr['PC2_1'] = PC2_1
	hdr['PC1_2'] = PC1_2
	hdr['PC2_2'] = PC2_2
	hdr['EQUINOX'] = 2000.0
	##write new fits
	hdul.writeto(new_fits_path+new_fits_filename+'.fits', overwrite=True)



def rdata(fits_path, fits_filename):
	with fits.open(fits_path+fits_filename+'.fits') as hdul:
		data0 = hdul[0].data
	return data0



def rdata3d(fits_path, fits_filename):
	with fits.open(fits_path+fits_filename+'.fits') as hdul:
		data0 = hdul[0].data
		data1 = hdul[1].data
	return data0, data1



def rakari(fits_path, fits_filename):
	with fits.open(fits_path+fits_filename+'.fits') as hdul:
		data0 = hdul[0].data
		data1 = hdul[1].data
		data2 = hdul[2].data
	return data0, data1, data2
