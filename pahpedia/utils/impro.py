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
	w = build_wcs(REFile)[0]
	data, footprint = reproject_interp(filIN+'.fits', hdr)
#	data, footprint = reproject_interp(read_fits(filIN, False), w, shape_out=(49,49)) # alternative
	comment = "Reprojected image."
	write_fits(filOUT, data, None, hdr, comment)#, \
#		BUNIT=hdr['BUNIT'], EQUINOX=hdr['EQUINOX'])

	return data, footprint, w

def crop(filIN, filOUT, centre, size, nopr=False):

	## read input.fits (2D ONLY)
	data, old_hdr = read_fits(filIN, False)
	hprint(nopr, "Raw size: ", data.shape)
	NAXIS1 = old_hdr['NAXIS1']
	NAXIS2 = old_hdr['NAXIS2']
	w, hdr = build_wcs(filIN)
	## convert coord
	w = WCS(old_hdr)
	xc, yc = w.all_world2pix(centre[0], centre[1], 1)
	if not (0<xc<NAXIS1 and 0<yc<NAXIS2):
		print("Error: crop center overpassed image border! ")
		exit()
	a = math.ceil(abs(size[0]/hdr['CDELT1'])) # ceil to avoid zero
	b = math.ceil(abs(size[1]/hdr['CDELT2']))
	hprint(nopr, "Crop size (input): ({}, {})".format(b,a))
	xmin = int(xc+math.floor(-a/2.)) # lowerleft centre
	xmax = int(xc+math.floor(a/2.))
	ymin = int(yc+math.floor(-b/2.))
	ymax = int(yc+math.floor(b/2.))
	hprint(nopr, "Crop size check (0 0 means OK): ", ymax-ymin-b, xmax-xmin-a)
	if not (xmin>=0 and xmax<=NAXIS1 and ymin>=0 and ymax<=NAXIS2):
		print("Error: crop region overpassed image border! ")
		exit()
	## modify header
	hdr['CRPIX1'], hdr['CRPIX2'] = a/2., b/2.
	hdr['CRVAL1'], hdr['CRVAL2'] = np.array(centre)
	## write output.fits
	comment = "Image cropped with centre ({:.5}, {:.5}) and size ({:.3}, {:.3}) (deg).".format(*centre, *size)
	write_fits(filOUT, data[ymin:ymax, xmin:xmax], None, hdr, comment)

	return data[ymin:ymax, xmin:xmax]

def cubislice(filIN, filOUT, *suffix):

	## read input.fits
	data, wvl, old_hdr = read_fits(filIN, True, 1)
	NAXIS3 = old_hdr['NAXIS3']
	for k in range(NAXIS3):
		## rebuild WCS object
		w, hdr = build_wcs(filIN, True)
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
	import matplotlib as mpl
	import matplotlib.pyplot as plt
	from myfunc import deg2world

	ref_path = '/Users/dhu/data/mosaic/SMC/'
	data_path = '../test_examples/'
	out_path = '../data/slices/'
	conv_path = '../data/convolved/'
	rpj_path = '../data/reprojection/'
	
	ref_filename = 'mips024'
	data_filename = 'n66_LL1_cube'
	out_ref = '_ref_'+ref_filename

	ra, dec = deg2world(0., 59., 3.5623, -72., 10., 33.972)
	d_ra, d_dec = 1./30., 1./30.

	ref = crop(ref_path+ref_filename, data_path+out_ref, \
		[ra, dec], [d_ra, d_dec])
	print("Cropped image size: ", ref.shape)
	## unit conversion
#	FLUXCONV = 0.000145730
#	ref = ref / FLUXCONV
	
	## 3D cube cropping
	## slice cube
	wvl = cubislice(data_path+data_filename, out_path+data_filename)
	filSL = []
	for k in range(np.size(wvl)):
		filSL.append(data_filename+'_'+'0'*(4-len(str(k)))+str(k))
	## crop slices
	cube=[]
	for sl_file in filSL:
		cube.append(crop(out_path+sl_file, out_path+sl_file, [ra, dec], [2.5*d_ra, 2.5*d_dec], nopr=True))
	cube = np.array(cube)
	print("Cropped cube size: ", cube.shape)

	## reprojection
	data, ft, w = rpj(conv_path+'n66_LL1_cube_0010_conv', rpj_path+'n66_LL1_cube_0010', data_path+out_ref)

	## plot cropped image
	plt.figure(figsize=(15,5))
	cmap = mpl.cm.viridis
	norm = mpl.colors.Normalize(vmin=0, vmax=1)
	ax1 = plt.subplot(131, projection=w)
	pcm = ax1.imshow(ref, origin='lower', cmap=cmap)

	ax2 = plt.subplot(132, projection=w)
	pcm = ax2.imshow(data, origin='lower', cmap=cmap)

	ax3 = plt.subplot(133, projection=w)
	pcm = ax3.imshow(ft, origin='lower', cmap=cmap)

	plt.colorbar(pcm, cmap=cmap)
	plt.show()
