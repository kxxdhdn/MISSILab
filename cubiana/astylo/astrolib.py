#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""

Astro library

"""
import sys, logging
logging.disable(sys.maxsize)

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS


def pc2cd(pc=None, cdelt=None, header=None, wcs=None):
	'''
	Convert CDELTia + PCi_ja to CDi_ja
	(astropy.wcs use PC/CDELT by default)

	------ INPUT ------
	pc                  PC matrix (priority if co-exist)
	cdelt               Coordinate increment at ref point
	header              header object (2nd priority if co-exist)
	wcs                 WCS object
	------ OUTPUT ------
	cl                  output object
	  cd                  CD matrix
	  pc                  PC matrix
	  cdelt               CDELTia
	'''
	## Initialize output object
	cl = type('', (), {})()
	cl.cd = np.zeros((2,2))
	cl.pc = np.zeros((2,2))
	cl.cdelt = np.zeros(2)

	if pc is not None and cdelt is not None:
		cl.pc = pc
		cl.cdelt = cdelt
		## CDi_j = PCi_j * CDELTi
		cl.cd = cl.pc * cl.cdelt.reshape((2,1))
	else:
		if header is not None:
			w = WCS(header)
		else:
			if WCS is not None:
				w = wcs
			else:
				print('ERROR: No input!')
				
				return cl

		cl.cd = w.pixel_scale_matrix

		if w.wcs.has_pc():
			cl.pc = w.wcs.get_pc()
			cl.cdelt = w.wcs.get_cdelt()

	return cl

# def cd2pc():
# 	cd = w.pixel_scale_matrix
# 	cdelt = np.sqrt((cd**2).sum(axis=0, dtype=float))
#	## See astropy.wcs.utils.proj_plane_pixel_scales

def fixwcs(file=None, header=None):
	'''
	Auto-detect & reduce dim if WCS is 3D with distortion

	------ INPUT ------
	file                FITS file (priority if co-exist)
	header              header object
	------ OUTPUT ------
	cl                  output object
	  header              header of primary HDU
	  wcs                 2D WCS
	  was3d               True: if input data is 3D
	'''
	## Initialize output object
	cl = type('', (), {})()
	cl.wcs = WCS(None, naxis=2)
	cl.header = None
	cl.was3d = False

	## Read file/header
	if file is not None:
		hdr = fits.open(file+'.fits')[0].header
		header = hdr.copy()
	else:
		if header is not None:
			hdr = header.copy()
		else:
			print('\nERROR: No input! \n')

			return cl

	## Reduce header dim/kw
	if header['NAXIS']==3:
		cl.was3d = True
		for kw in hdr.keys():
			if '3' in kw:
				del header[kw]
		header['NAXIS'] = 2
		header['COMMENT'] = "3D keywords excluded (for astropy.wcs). "
	
	## Create 2D WCS object
	cl.wcs = WCS(header, naxis=2)

	cl.header = header # (reduced) header

	return cl

