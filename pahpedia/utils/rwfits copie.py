#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""

Read & write .fits file

"""

from astropy.io import fits
import numpy as np

path = "./"
#path = "/Users/abricot/Github/"

def read_fits(filename, hdr, data_on=True, wvl_on=True, w_mod=False):
	with fits.open(path+filename+'.fits') as hdul:
		if data_on==True:
			data = hdul[0].data

		if wvl_on==True:
			if w_mod==0:
				wvl = hdul[1].data[0][0][:,0]
			else:
				wvl = hdul[1].data # for rewitten header

		value = []
		for name in range(np.size(hdr)):
			value.append(hdul[0].header[hdr[name]])
			

	if data_on==True:
		if wvl_on==True:
			return data, wvl, value
		else:
			return data, value
	else:
		if wvl_on==True:
			return wvl, value
		else:
			return value

#def write_fits():


"""
------------------------------ MAIN (test) ------------------------------
"""
if __name__ == "__main__":

	## in situ test
	path = '../test_examples/'
	filename = 'n66_LL1_cube'

	data, wvl, hdr = read_fits(path+filename, ['NAXIS3'])

	print(data, wvl, hdr)
