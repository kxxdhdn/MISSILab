#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS

global fitsext
fitsext = '.fits'

def read(file, file_unc=None, wmod=0):
	'''
	Read fits file (auto detect dim)

	------ INPUT ------
	file                input fits file
	file_unc            input uncertainty file
	wmod                output wave mode (Default: 0)
	                      0 - 1darray; 
	                      1 - FITS_rec.
	------ OUTPUT ------
	dataset             dataset object
	  header              header of primary HDU
	  data                data in primary HDU
	  header_w            header of W-TAB
	  wave                data in table 1 (None if does not exist)
	  unc                 uncertainty array
	'''
	## Initialize dataset object
	dataset = type('', (), {})()
	dataset.header_w = None
	dataset.wave = None
	dataset.unc = None

	## Read header & data
	with fits.open(file+fitsext) as hdul:
		hdr = hdul[0].header
		dataset.data = hdul[0].data
		dataset.header = hdr
		if 'CTYPE3' in hdr.keys():
			del hdr['CTYPE3']
		dataset.WCS = WCS(hdr)

		## Read wavelength
		if len(hdul)==2:
			dataset.header_w = hdul[1].header
			wave = hdul[1].data

			if isinstance(hdul[1], fits.BinTableHDU):
				if wmod==0:
					wave = wave[0][0][:,0] ## Convert FITS_rec to 1darray
			elif isinstance(hdul[1], fits.ImageHDU):
				Nw = len(wave)
				if wmod==1:
					wave = np.array(wave).reshape((Nw,1))
					col = fits.Column(array=[wave], format=str(Nw)+'E', \
						name='WAVE-TAB', unit='um', dim='(1,{})'.format(Nw))
					tab = fits.BinTableHDU.from_columns([col], name='WCS-TAB ')
					wave = tab.data
			
			dataset.wave = wave
	
	if file_unc is not None:
		## Read uncertainty data
		with fits.open(file+fitsext) as hdul:
			dataset.unc = hdul[0].data
		
	return dataset

def red_wcs(file=None, header=None):
	'''
	Extract WCS (auto detect & reduce dim if 3D)

	------ INPUT ------
	file                input fits file (priority if co-exist with input header)
	header              input fits header
	------ OUTPUT ------
	dataset             dataset object
	  header              header of primary HDU
	  WCS                 2D WCS
	  is3d                True: if input data is 3D
	'''
	## Initialize dataset object
	dataset = type('', (), {})()
	dataset.WCS = WCS(None, naxis=2)
	dataset.header = None
	dataset.is3d = False

	## Read file/header
	if file is not None:
		hdr = fits.open(file+fitsext)[0].header
		header = hdr.copy()
	else:
		if header is not None:
			hdr = header.copy()
		else:
			return dataset

	## Reduce header dim/kw
	if header['NAXIS']==3:
		dataset.is3d = True
		for kw in hdr.keys():
			if '3' in kw:
				del header[kw]
		header['NAXIS'] = 2
		header['COMMENT'] = "3D distorsion kw excluded (for astropy.wcs). "
	
	## Create 2D WCS object
	dataset.WCS = WCS(header, naxis=2)

	dataset.header = header # (reduced) header

	return dataset

def write_fits(file, header, data, wave=None, wmod=0, **hdrl):
	'''
	Write fits file

	------ INPUT ------
	file                input fits file
	header              header of primary HDU
	data                data in primary HDU
	wave                data in table 1 (ndarray. default: None)
	wmod                wave table format (0 - Image; 1 - BinTable. default: 0)
	------ OUTPUT ------
	'''
	for key, value in hdrl.items():
		header[key] = value
	primary_hdu = fits.PrimaryHDU(header=header, data=data)
	hdul = fits.HDUList(primary_hdu)
	
	## Add table
	if wave is not None:
		## Convert wave format
		if isinstance(wave, fits.fitsrec.FITS_rec):
			if wmod==0:
				wave = wave[0][0][:,0]
		else:
			Nw = len(wave)
			if wmod==1:
				wave = np.array(wave).reshape((Nw,1))
				col = fits.Column(array=[wave], format=str(Nw)+'E', \
					name='WAVE-TAB', unit='um', dim='(1,{})'.format(Nw))
				tab = fits.BinTableHDU.from_columns([col], name='WCS-TAB ')
				wave = tab.data
		## Create table
		if wmod==0:
			hdu = fits.ImageHDU(data=wave, name='WAVE-TAB')
		elif wmod==1:
			hdu = fits.BinTableHDU(data=wave, name='WCS-TAB ')

		hdul.append(hdu)

	hdul.writeto(file+fitsext, overwrite=True)
	