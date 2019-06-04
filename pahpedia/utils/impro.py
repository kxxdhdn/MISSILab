#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""

PSF homogeneisation

"""

from astropy.io import ascii
import numpy as np
import math
from reproject import reproject_interp
import subprocess as SP
from rwfits import *

mu, sigma = 0, 1.

def specorrect(filIN, filOUT, factor, wmin, wmax, wmod=0):
	
	im, wvl, hdr = read_fits(filIN)
	for k, lam in enumerate(wvl):
		if lam>wmin and lam<wmax:
			im[k,:,:] *= factor
	write_fits(filOUT, im, hdr, wvl)

	return im

def swarpj(filIN, filOUT, REFile):

	w, hdr = WCSextract(REFile)
	CRVAL1 = str(hdr['CRVAL1'])
	CRVAL2 = str(hdr['CRVAL2'])
	f = open('coadd.head', 'w')
	for i, card in enumerate(hdr):
		if card!='':
			f.write(card+' = '+str(hdr[i])+'\n')
	f.close()

	SP.call('swarp -d >swarp.cfg', shell=True)
	SP.call('swarp '+filIN+'.fits -c swarp.cfg \
		-IMAGEOUT_NAME '+filOUT+'.fits \
		-RESAMPLING_TYPE BILINEAR -SUBTRACT_BACK N \
		-CENTER_TYPE MANUAL,MANUAL -CENTER '+CRVAL1+','+CRVAL2,
		shell=True)

	data = read_fits(filOUT, False)[0]

	return data, w

def rpj(filIN, filOUT, hdREF, write=True):

	data, footprint = reproject_interp(filIN+'.fits', hdREF)
	comment = "Reprojected image."
	if write==True:
		write_fits(filOUT, data, hdREF, COMMENT=comment)

	return data, footprint

def hextract(filIN, filOUT, x0, x1, y0, y1):
	
	oldim = read_fits(filIN, False)[0]
	w, hdr = WCSextract(filIN)
	# hdr['NAXIS1'] = x1 - x0 + 1
	# hdr['NAXIS2'] = y1 - y0 + 1
	hdr['CRPIX1'] += -x0
	hdr['CRPIX2'] += -y0
	newim = oldim[y0:y1+1, x0:x1+1]

	write_fits(filOUT, newim, hdr)

	return newim

def crop(filIN, filOUT, cen, size):

	print("Crop centre (ra, dec): [{:.8}, {:.8}]".format(*cen))
	print("Crop size (pix): [{}, {}].".format(*size))

	## read input.fits (2D ONLY)
	data = read_fits(filIN, False)[0]
	print("Raw size (pixel): ", data.shape)
	w, hdr = WCSextract(filIN)
	NAXIS1 = hdr['NAXIS1']
	NAXIS2 = hdr['NAXIS2']

	## convert coord
	xc, yc = w.all_world2pix(cen[0], cen[1], 1)
	if not (0<xc<NAXIS1 and 0<yc<NAXIS2):
		print("Error: crop centre overpassed image border! ")
		exit()
	## lowerleft origin
	xl = math.floor(xc - size[0]/2.)
	yl = math.floor(yc - size[1]/2.)
	xr = xl + size[0]
	yr = yl + size[1]
	if not (xl>=0 and xr<=NAXIS1 and yl>=0 and yr<=NAXIS2):
		print("Error: crop region overpassed image border! ")
		exit()
	## modify header
	hdr['CRPIX1'] += -xl
	hdr['CRPIX2'] += -yl
	## write output.fits
	comment = "Image cropped at centre: [{:.8}, {:.8}]. ".format(*cen)
	comment = "with size [{}, {}] (pixels).".format(*size)
	write_fits(filOUT, data[yl:yr, xl:xr], hdr, COMMENT=comment)

	return data[yl:yr, xl:xr]

def cubislice(filIN, filOUT, uncIN=None, filOFF=None, wmod=0, *suffix):

	## read input.fits
	im, wvl = read_fits(filIN, True, wmod)[0:2]
	hdr = WCSextract(filIN)[1]

	if uncIN!=None:
		unc = read_fits(uncIN)[0]
		NAXIS1 = hdr['NAXIS1']
		NAXIS2 = hdr['NAXIS2']

	if filOFF!=None:
		back = read_fits(filOFF, True, wmod)[0]
		off = np.nanmean(back, axis=(1,2))

	for k in range(np.size(wvl)):

		if uncIN!=None:
			theta = np.random.normal(mu, sigma, NAXIS1*NAXIS2).reshape(NAXIS2, NAXIS1)
			im[k,:,:] = im[k,:,:] + theta * unc[k,:,:]
		## output filename list
		filSL = filOUT+'_'+'0'*(4-len(str(k)))+str(k)

		if filOFF!=None:
			im[k,:,:] += -off[k]

		for s in suffix:
			filSL = filSL+s
		
		# comment = "NO.{} image sliced from {}".format(k, filIN)
		write_fits(filSL, im[k,:,:], hdr)

	return wvl, hdr
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
	out_path = '../test_examples/slices/'
	
	ref_filename = 'mips024'
	data_filename = 'n66_LL1'
	out_ref = '_ref_'+ref_filename

	ra, dec = celest2deg(0., 59., 3.5623, -72., 10., 33.972)
	dx, dy = 34, 40

	ref = crop(ref_path+ref_filename, data_path+out_ref, \
		(ra, dec), (dx, dy))
	print("Cropped image size: ", ref.shape)

	hdr = read_fits(data_path+out_ref, False)[1]
	
	## 3D cube cropping
	## slice cube
	wvl = cubislice(data_path+data_filename, out_path+data_filename, None, None, 1)
	## reprojection
	data_nocrop = rpj(out_path+'n66_LL1_0010', data_path+'n66_LL1_cube_0010_nocrop', hdr)[0]
	exit()
	
	filSL = []
	for k in range(np.size(wvl)):
		filSL.append(data_filename+'_'+'0'*(4-len(str(k)))+str(k))
	## crop slices
	cube=[]
	for sl_file in filSL:
		cube.append(crop(out_path+sl_file, out_path+sl_file, (dec, ra), (dy, dx)))
	cube = np.array(cube)
	print("Cropped cube size: ", cube.shape)

	## reprojection
	data, ft = rpj(out_path+'n66_LL1_0010', data_path+'n66_LL1_0010', hdr)
	
	## crop match test
	imview([ref, data, data-ref, data-data_nocrop], figshape=(1,4), figsize=(15,5))
	# imview([data, data], w, (1,2), figsize=(12,5))
	# imview([data - data_nocrop, data, data_nocrop], w, (1,3), figsize=(12,5))

	plt.show()
