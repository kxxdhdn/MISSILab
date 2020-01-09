#!/usr/bin/env python
# -*- coding: utf-8 -*-

import time
t0 = time.time()

import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))+'/..')

from astropy.io import fits
from astropy.wcs import WCS

from astylo.bio import read_fits, ext_wcs

## TEST read_fits
##----------------
# FITSdata = read_fits('data/M83_LL1')
# print(FITSdata.WCS)

## TEST WCS
##----------
hdul = fits.open('data/M83_LL1.fits')
hdr = hdul[0].header
header = hdr.copy()
# for kw in hdr.keys():
	# if '3' in kw:
		# del header[kw]
		# print(kw)

# del header['NAXIS3']
# del header['PC3_3']
# del header['CRPIX3']
# del header['CRVAL3']
del header['CTYPE3'] # the one (with which SIP kw counted)
# del header['CUNIT3']
# del header['PS3_0']
# del header['PS3_1']
# del header['FLXCON03']
# del header['FLXERR03']
# del header['GAIN3']
# del header['BASECH3']

# print(WCS(header))
# print(list(header.keys()))

# print(WCS('data/M83_LL1.fits')) # not working

# w = ext_wcs('data/M83_LL1').WCS
w = ext_wcs().WCS
print(w[1])

t_total = time.time()
print("****** total_time = {:.0f} seconds ******".format(t_total - t0))
