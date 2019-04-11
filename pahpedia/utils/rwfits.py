#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""

Read & write .fits file

"""

from astropy.io import fits
from astropy.wcs import WCS
import numpy as np

path = "./"
#path = "/Users/abricot/Github/"

def read_fits(filename, hdrl, data_on=True, wvl_on=True, w_mod=False):

	with fits.open(path+filename+'.fits') as hdul:
		if data_on==True:
			data = hdul[0].data
		else:
			data = None

		if wvl_on==True:
			if w_mod==0:
				wvl = hdul[1].data[0][0][:,0]
			else:
				wvl = hdul[1].data # for rewitten header
		else:
			wvl = None

		hdr = hdul[0].header
		#wd = WCS(hdr)
		value = []
		for i in range(np.size(hdrl)):
			value.append(hdr[hdrl[i]])
			
	return data, wvl, value#, wd

def write_fits(filename, data, wvl, comment, **hdrl):

	hdr = fits.Header()
	for key, value in hdrl.items():
		hdr[key] = value
	hdr['EQUINOX'] = 2000.0
	primary_hdu = fits.PrimaryHDU(header=hdr, data=data)
	hdul = fits.HDUList(primary_hdu)
	## add table
	hdu = fits.ImageHDU(data=wvl, name="Wavelength (microns)")
	hdul.append(hdu)

	hdul.writeto(path+filename+'.fits', overwrite=True)



"""
------------------------------ MAIN (test) ------------------------------
"""
if __name__ == "__main__":

	## in situ test
	path = '../test_examples/'
	filename = 'n66_LL1_cube'
	writename = 'write_test'

	data, wvl, hdrl = read_fits(path+filename, ['NAXIS3', 'NAXIS2', 'NAXIS1'])
	#print(data, wvl, hdrl)

	write_fits(path+writename, data, wvl, "This is a write test", 
		NAXIS1=hdrl[2], NAXIS2=hdrl[1], NAXIS3=hdrl[0])
	data_, wvl_, hdrl_ = read_fits(path+writename, ['NAXIS3', 'NAXIS2', 'NAXIS1'], w_mod=True)
	print(data_, wvl_, hdrl_)
