#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""

Read & write .fits file

"""

from astropy.io import fits
from astropy.wcs import WCS
import numpy as np

def read_fits(filename, is3d=True, wvl_mod=0, data_on=True):

	with fits.open(filename+'.fits') as hdul:
		## read header
		hdr = hdul[0].header

		## read data
		if is3d==True: # 3D
			if wvl_mod==0:
				wvl = hdul[1].data # rewitten header
			else:
				wvl = hdul[1].data[0][0][:,0] # CUBISM witten
			if data_on==True:
				data = hdul[0].data
				return data, wvl, hdr
			else:
				return wvl, hdr
		else: # 2D
			if data_on==True:
				data = hdul[0].data
				return data, hdr
			else:
				return hdr

def WCSextract(filename):

	hdr = fits.open(filename+'.fits')[0].header

	hdr0 = hdr.copy()
	if hdr['NAXIS']==3:
		for kw in hdr0.keys():
			if '3' in kw:
				del hdr[kw]
		hdr['COMMENT'] = "This header is adapted to 2D WCS extraction need. "
	
	w = WCS(hdr, naxis=2)

	return w, hdr

def write_fits(filename, data, wvl, hdr, **hdrl):

	for key, value in hdrl.items():
		hdr[key] = value
	primary_hdu = fits.PrimaryHDU(header=hdr, data=data)
	hdul = fits.HDUList(primary_hdu)
	## add table
	if wvl is not None:
		hdu = fits.ImageHDU(data=wvl, name="Wavelength (microns)")
		hdul.append(hdu)

	hdul.writeto(filename+'.fits', overwrite=True)

"""
------------------------------ MAIN (test) ------------------------------
"""
if __name__ == "__main__":

	## in situ test
	path = '../test_examples/'
	filename = 'n66_LL1_cube'
	writename = 'write_test'

	data, wvl, hdr = read_fits(path+filename, True, 1)
	"""
	print(data, wvl, hdr)
	w, new_hdr = WCSextract(path+filename)
	print(w)
	write_fits(path+writename, data, wvl, hdr, \
		COMMENT="This is a write test", NAXIS3=hdr['NAXIS3'])
	data_, wvl_, new_hdr = read_fits(path+writename)
	print(data_, wvl_, new_hdr)
	"""
	## 2D write
	w, new_hdr = WCSextract(path+filename)
	# print(w)
	write_fits(path+writename, data[0,:,:], None, new_hdr, \
		COMMENT="This is a 2D write test", \
		BUNIT=hdr['BUNIT'], EQUINOX=hdr['EQUINOX'])
	data_, new_hdr = read_fits(path+writename, False)
	# print(data_, new_hdr)
	# print(new_hdr)
