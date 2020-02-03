#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys, os, logging
testdir = os.path.dirname(os.path.abspath(__file__))
logging.disable(sys.maxsize)

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS

## Local
sys.path.insert(0, testdir+'/..') ## astylo path
from astylo.astrolib import fixwcs, pc2cd

## Set path
datdir = testdir+'/dat/'
outdir = testdir+'/out/'


## TEST pc2cd
##------------
pc = np.array([[-0.695315298625, -0.71870483197], [0.71870483197, -0.695315298625]])
cdelt = np.array([-0.000513888895512, 0.000513888895512])
print(pc2cd(pc, cdelt).cd, '\n')

hdr = fixwcs(datdir+'M82_09_SL2').header
print(pc2cd(header=hdr).cd)
print(pc2cd(header=hdr).pc)
print(pc2cd(header=hdr).cdelt, '\n')

w = fixwcs(datdir+'M82_09_SL2').wcs
print(pc2cd(wcs=w).cd, '\n')

'''

## TEST fixwcs
##-------------
# print(WCS(datdir+'M82_09_LL2.fits')) # not working

# hdul = fits.open('data/M83_LL1.fits')
hdr = fits.open(datdir+'M82_09_SL2.fits')[0].header
header = hdr.copy()
# for kw in hdr.keys():
	# if '3' in kw:
		# del header[kw]
		# print(kw)

header['NAXIS'] = 2 # SIP kw sensible
del header['NAXIS3']
del header['PC3_3'] # SIP kw sensible
del header['CRPIX3'] # SIP kw sensible
del header['CRVAL3'] # SIP kw sensible
del header['CTYPE3'] # SIP kw sensible
del header['CUNIT3'] # SIP kw sensible
del header['PS3_0'] # SIP kw sensible
del header['PS3_1'] # SIP kw sensible
# del header['FLXCON03']
# del header['FLXERR03']
# del header['GAIN3']
# del header['BASECH3']
print(WCS(header))
# print(list(header.keys()))

# print(fixwcs(datdir+'M82_09_LL2').header)
print(fixwcs().wcs)

'''