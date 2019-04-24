#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""

PSF homogeneisation

"""

import numpy as np
import math
from reproject import reproject_interp
from rwfits import *
from myfunc import hprint

def rpj(filIN, filOUT, REFile):

	## 2D ONLY
	hdr = read_fits(REFile, False)[1]
	w = WCSextract(REFile)[0]
	data, footprint = reproject_interp(filIN+'.fits', hdr)
#	data, footprint = reproject_interp(read_fits(filIN, False), w, shape_out=(49,49)) # alternative
	comment = "Reprojected image."
	write_fits(filOUT, data, None, hdr, comment)#, \
#		BUNIT=hdr['BUNIT'], EQUINOX=hdr['EQUINOX'])

	return data, footprint, w

def crop(filIN, filOUT, CROrigin, CROsize, nopr=False):

	## read input.fits (2D ONLY)
	data, old_hdr = read_fits(filIN, False)
	NAXIS1 = old_hdr['NAXIS1']
	NAXIS2 = old_hdr['NAXIS2']
	hprint(nopr, "Raw size: ", (NAXIS2, NAXIS1))

	w, hdr = WCSextract(filIN)
	## convert coord
	x0, y0 = w.wcs_world2pix(CROrigin[0], CROrigin[1], 1)
	x0 = math.floor(x0)
	y0 = math.floor(y0)
	hprint(nopr, "y0, x0 = ", y0, x0)
	if not (0<=x0<NAXIS1 and 0<=y0<NAXIS2):
		print("Error: crop lowerleft overpassed image border! ")
		exit()
	a = math.ceil(abs(CROsize[0]/hdr['CDELT1'])) # ceil to avoid zero
	b = math.ceil(abs(CROsize[1]/hdr['CDELT2']))
	hprint(nopr, "Crop size (input): ({}, {})".format(b,a))
	xmin = x0
	xmax = x0 + a
	ymin = y0
	ymax = y0 + b
	if not (xmax<=NAXIS1 and ymax<=NAXIS2):
		print("Error: crop region overpassed image border! ")
		exit()
	## modify header
	hdr['CRPIX1'], hdr['CRPIX2'] = a/2., b/2.
	hdr['CRVAL1'], hdr['CRVAL2'] = np.array(CROrigin) + np.array(CROsize) / 2.
	## write output.fits
	comment = "Image cropped with lowerleft origin ({:.5}, {:.5}) and size ({:.3}, {:.3}) (deg).".format(*CROrigin, *CROsize)
	write_fits(filOUT, data[ymin:ymax, xmin:xmax], None, hdr, comment)

	return data[ymin:ymax, xmin:xmax]

def cubislice(filIN, filOUT, *suffix):

	## read input.fits
	data, wvl, old_hdr = read_fits(filIN, True, 1)
	NAXIS3 = old_hdr['NAXIS3']
	for k in range(NAXIS3):
		## rebuild WCS object
		w, hdr = WCSextract(filIN)
		filSL = filOUT+'_'+'0'*(4-len(str(k)))+str(k)
		for s in suffix:
			filSL = filSL+s
		comment = "NO.{} image sliced from {}".format(k, filIN)
		write_fits(filSL, data[k,:,:], wvl, hdr, comment, \
			BUNIT=old_hdr['BUNIT'], EQUINOX=old_hdr['EQUINOX'])

	return wvl
"""
------------------------------ MAIN (test) ------------------------------
"""
if __name__ == "__main__":

	## in situ test
	from astropy import units as u
	import matplotlib.pyplot as plt
	from myfunc import celest2deg
	from mapping import *

	ref_path = '/Users/dhu/data/mosaic/SMC/'
	data_path = '../test_examples/'
	out_path = '../data/slices/'
	conv_path = '../data/convolved/'
	rpj_path = '../data/reprojection/'
	
	ref_filename = 'mips024'
	data_filename = 'n66_LL1_cube'
	out_ref = '_ref_'+ref_filename

	ra, dec = celest2deg(0., 59., 3.5623, -72., 10., 33.972)
	d_ra, d_dec = 1./30., 1./25.
	print("d_ra, d_dec", d_ra, d_dec)

	ref = crop(ref_path+ref_filename, data_path+out_ref, \
		(ra, dec), (d_ra, d_dec), nopr=0)
	print("Cropped image size: ", ref.shape)
	## unit conversion
	FLUXCONV = 0.000145730
	ref = ref / FLUXCONV
	data, hdr = read_fits(data_path+out_ref, False)
	hdr['PC1_1'], hdr['PC1_2'] = 1., 0.
	hdr['PC2_1'], hdr['PC2_2'] = 0., 1.
	comment="add pc values"
	write_fits(data_path+out_ref, data, None, hdr, comment)
	
	## 3D cube cropping
	## slice cube
	wvl = cubislice(data_path+data_filename, out_path+data_filename)
	data_nocrop = rpj(out_path+'n66_LL1_cube_0010', rpj_path+'n66_LL1_cube_0010_nocrop', data_path+out_ref)[0]

	filSL = []
	for k in range(np.size(wvl)):
		filSL.append(data_filename+'_'+'0'*(4-len(str(k)))+str(k))
	## crop slices
	cube=[]
	for sl_file in filSL:
		cube.append(crop(out_path+sl_file, out_path+sl_file, [ra, dec], [1.*d_ra, 1.*d_dec], nopr=1))
	cube = np.array(cube)
	print("Cropped cube size: ", cube.shape)

	## reprojection
	data, ft, w = rpj(out_path+'n66_LL1_cube_0010', rpj_path+'n66_LL1_cube_0010', data_path+out_ref)
	ref = rpj(data_path+out_ref, rpj_path+out_ref, data_path+out_ref)[0]
	## plot cropped image
	# imview([ref, data, ft], w, (1,3), figsize=(12,5))

	## crop match test
	imview([ref, data, data_nocrop], w, (1,3), figsize=(12,5))

	plt.show()
