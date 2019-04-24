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

def rpj(filIN, filOUT, REFile, size):

	## 2D ONLY
	w, hdr = WCSextract(REFile)
	data, footprint = reproject_interp(filIN+'.fits', w, shape_out=size)
	comment = "Reprojected image."
	write_fits(filOUT, data, None, hdr, COMMENT=comment)

	return data, footprint, w

def crop(filIN, filOUT, centre, size, nopr=False):

	## read input.fits (2D ONLY)
	data = read_fits(filIN, False)[0]
	hprint(nopr, "Raw size: ", data.shape)
	w, hdr = WCSextract(filIN)
	NAXIS1 = hdr['NAXIS1']
	NAXIS2 = hdr['NAXIS2']
	## convert coord
	xc, yc = w.all_world2pix(centre[1], centre[0], 1)
	xc0 = math.floor(xc)
	yc0 = math.floor(yc)
	ra0, dec0 = w.all_pix2world(xc0, yc0, 1)
	# print("---------------", ra0-centre[0], dec0-centre[1])
	hprint(nopr, "Crop centre (input): [{:.5}, {:.5}]".format(*centre))
	hprint(nopr, "Crop centre: [{:.5}, {:.5}]".format(dec0, ra0))
	if not (0<xc<NAXIS1 and 0<yc<NAXIS2):
		print("Error: crop center overpassed image border! ")
		exit()
	hprint(nopr, "Crop size (input): ", size)
	## lowerleft centre
	xmin = xc0 - math.floor(size[1]/2.)
	xmax = xc0 + math.ceil(size[1]/2.)
	ymin = yc0 - math.floor(size[0]/2.)
	ymax = yc0 + math.ceil(size[0]/2.)
	if not (xmin>=0 and xmax<=NAXIS1 and ymin>=0 and ymax<=NAXIS2):
		print("Error: crop region overpassed image border! ")
		exit()
	## modify header
	hdr['CRPIX1'], hdr['CRPIX2'] = size[1]/2., size[0]/2.
	hdr['CRVAL1'], hdr['CRVAL2'] = float(ra0), float(dec0)
	## write output.fits
	comment = "Image cropped at centre ({:.5}, {:.5}) with pixel size ({}, {}).".format(ra0, dec0, size[1], size[0])
	write_fits(filOUT, data[ymin:ymax, xmin:xmax], None, hdr, COMMENT=comment)

	return data[ymin:ymax, xmin:xmax]

def cubislice(filIN, filOUT, *suffix):

	## read input.fits
	data, wvl = read_fits(filIN, True, 1)[0:2]
	hdr = WCSextract(filIN)[1]
	for k in range(np.size(wvl)):
		filSL = filOUT+'_'+'0'*(4-len(str(k)))+str(k)
		for s in suffix:
			filSL = filSL+s
		
		# comment = "NO.{} image sliced from {}".format(k, filIN)
		write_fits(filSL, data[k,:,:], wvl, hdr)

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
	dx, dy = 34, 40

	ref = crop(ref_path+ref_filename, data_path+out_ref, \
		(dec, ra), (dy, dx))
	print("Cropped image size: ", ref.shape)
	## unit conversion
	# FLUXCONV = 0.000145730
	# ref = ref / FLUXCONV
	
	## 3D cube cropping
	## slice cube
	wvl = cubislice(data_path+data_filename, out_path+data_filename)
	## reprojection
	data_nocrop = rpj(out_path+'n66_LL1_cube_0010', rpj_path+'n66_LL1_cube_0010_nocrop', data_path+out_ref, (dy, dx))[0]

	filSL = []
	for k in range(np.size(wvl)):
		filSL.append(data_filename+'_'+'0'*(4-len(str(k)))+str(k))
	## crop slices
	cube=[]
	for sl_file in filSL:
		cube.append(crop(out_path+sl_file, out_path+sl_file, (dec, ra), (dy, dx), nopr=True))
	cube = np.array(cube)
	print("Cropped cube size: ", cube.shape)

	## reprojection
	data, ft, w = rpj(out_path+'n66_LL1_cube_0010', rpj_path+'n66_LL1_cube_0010', data_path+out_ref, (dy, dx))
	# print("w", w)
	# print("wref", WCSextract(ref_path+ref_filename)[0])
	# print("w0", WCSextract(data_path+data_filename)[0])
	# print("w10", WCSextract(out_path+'n66_LL1_cube_0010')[0])
	
	## crop match test
	multimview([ref, data, data-ref, data-data_nocrop], w, (1,4), figsize=(15,5))
	# imview([data, data], w, (1,2), figsize=(12,5))
	# imview([data - data_nocrop, data, data_nocrop], w, (1,3), figsize=(12,5))

	plt.show()
