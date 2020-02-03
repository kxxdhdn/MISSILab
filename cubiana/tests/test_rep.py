#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))+'/..')

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from reproject import reproject_interp
import matplotlib.pyplot as plt

from astylo.bio import write_fits


## imontage rep
##----------------------------------------------
# hdul1 = fits.open('data/M82_09_SL1_rep.fits')
# hdr1 = hdul1[0].header
# header = hdr1.copy()
# header['NAXIS'] = 2 #
# for kw in hdr1.keys():
# 	if '3' in kw:
# 		del header[kw]
# header['CTYPE1'] = 'RA---TAN-SIP'
# header['CTYPE2'] = 'DEC--TAN-SIP'
# w1 =WCS(header)
# im_rep = hdul1[0].data[0]


## REF
##----------------------------------------------
# hdul0 = fits.open('data/footprint0.fits')
hdul0 = fits.open('data/M82_SL2_fp.fits')
hdr0 = hdul0[0].header
header = hdr0.copy()
header['NAXIS'] = 2 #
for kw in hdr0.keys():
	if '3' in kw:
		del header[kw]
header['CTYPE1'] = 'RA---TAN-SIP'
header['CTYPE2'] = 'DEC--TAN-SIP'
w0 = WCS(header)
hdr_out = w0.to_header(relax=True)


## Origin
##----------------------------------------------
hdu = fits.open('data/M82_09_SL1.fits')[0]

hdr = hdu.header
header = hdr.copy()
header['NAXIS'] = 2 #
# del header['NAXIS3']
del header['PC3_3'] #
del header['CRPIX3'] #
del header['CRVAL3'] #
del header['CTYPE3'] #
del header['CUNIT3'] #
del header['PS3_0'] #
del header['PS3_1'] #
# for kw in hdr.keys():
# 	if '3' in kw:
# 		del header[kw]
header['CTYPE1'] = 'RA---TAN-SIP'
header['CTYPE2'] = 'DEC--TAN-SIP'
w = WCS(header)
w = w.sub((1,2))
hdu.header = w.to_header(relax=True)

s0 = hdu.data[0]
hdu.data = s0
shape_out = hdr0['NAXIS2'], hdr0['NAXIS1']

# im = reproject_interp(hdu, w0, 
# 	shape_out=shape_out)[0]
im = reproject_interp((s0, w), w0, 
	shape_out=shape_out)[0]

write_fits('./out/test_rep',hdr_out, im)


## Imshow check
##----------------------------------------------
pix = [[32,12], [32,48], [64,48], [96,12], [96,48], \
	[96,84], [124,48], [160,48], [160,84]]
for p in pix:
	s0[p[1], p[0]] = np.nan
world = w.all_pix2world(np.array(pix), 1)
pix_rep = w0.all_world2pix(np.array(world), 1)
# print(pix_rep)
for p in pix_rep:
	x = int(p[0])
	y = int(p[1])
	# print(p)
	im[y, x] = np.nan

# p0 = [85,70]
# s0[p0[1]-5:p0[1]+6, p0[0]-6:p0[0]+6] = np.nan
# ra, dec = w.all_pix2world(p0[0], p0[1], 1)
# x, y = w1.all_world2pix(ra, dec, 1)
# x = int(x)
# y = int(y)
# im_rep[y-5:y+6, x-5:x+6] = np.nan

x0 = np.arange(80, 91)
y0 = np.arange(65, 76)
s0[y0[0]:y0[-1], x0[0]:x0[-1]] = np.nan
p0 = []
for i in range(10):
	for j in range(10):
		p0.append([x0[i], y0[j]])
world0 = w.all_pix2world(np.array(p0), 1)
p0_rep = w0.all_world2pix(world0, 1)
for p in p0_rep:
	x = int(p[0])
	y = int(p[1])
	im[y, x] = np.nan


## Plot
##----------------------------------------------
fig=plt.figure()
col = 2
row = 1
# for i in range(1, col*row+1):
## Test
##------
fig.add_subplot(row, col, 1, projection=w0)
plt.imshow(im, cmap=plt.cm.viridis)
plt.colorbar()
# for p in pix_rep:
# 	plt.scatter(p[0], p[1], color='r')
## Origin
##--------
fig.add_subplot(row, col, 2, projection=w)
plt.imshow(s0, cmap=plt.cm.viridis)
plt.colorbar()
# for p in pix:
# 	plt.scatter(p[0], p[1], color='r')
## Reprojection
##--------------
# fig.add_subplot(row, col, 2, projection=w0)
# plt.imshow(im, cmap=plt.cm.viridis)
# plt.colorbar()
# for p in pix_rep:
# 	plt.scatter(p[0], p[1], color='r')

plt.show()
