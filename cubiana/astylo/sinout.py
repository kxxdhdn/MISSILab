#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""

Sublime input & output

"""
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
import h5py as H5
import csv

global fitsext, h5ext
fitsext = '.fits'
h5ext = '.h5'
ascext = '.txt'
csvext = '.csv'

def read_fits(file, wmod=0):
	'''
	Read fits file (auto detect dim)

	--- INPUT ---
	file        input fits file
	wmod        wave mode
	--- OUTPUT ---
	hdr         header of primary HDU
	data        data in primary HDU
	wave        data in table 1 (if exists)
	'''
	with fits.open(file+fitsext) as hdul:
		## read header
		hdr = hdul[0].header

		## read data
		# 3D
		if hdr['NAXIS']==3:
			data = hdul[0].data
			
			if wmod==0:
				wave = hdul[1].data # rewitten header
			else:
				wave = hdul[1].data[0][0][:,0] # CUBISM witten
			
			return hdr, data, wave
		# 2D
		else:
			data = hdul[0].data

			return hdr, data

def write_fits(file, hdr, data, wave=None, **hdrl):
	'''
	Write fits file

	--- INPUT ---
	file        input fits file
	hdr         header of primary HDU
	data        data in primary HDU
	wave        data in table 1 (default: None)
	--- OUTPUT ---
	new fits file
	'''
	for key, value in hdrl.items():
		hdr[key] = value
	primary_hdu = fits.PrimaryHDU(header=hdr, data=data)
	hdul = fits.HDUList(primary_hdu)
	## add table
	if wave is not None:
		hdu = fits.ImageHDU(data=wave, name="Wavelength (microns)")
		hdul.append(hdu)

	hdul.writeto(file+fitsext, overwrite=True)

def WCSextract(file):
	'''
	extract WCS (auto detect & reduce dim if 3D)

	--- INPUT ---
	file        input fits file
	--- OUTPUT ---
	hdr         header of primary HDU
	w           2D WCS
	is3d        if input data is 3D: True
	'''
	hdr0 = fits.open(file+fitsext)[0].header

	hdr = hdr0.copy()
	if hdr['NAXIS']==3:
		is3d = True
		for kw in hdr0.keys():
			if '3' in kw:
				del hdr[kw]
		hdr['NAXIS'] = 2
		hdr['COMMENT'] = "This header is adapted to 2D WCS extraction need. "
	else:
		is3d = False
	
	w = WCS(hdr, naxis=2)

	return hdr, w, is3d

def read_hdf5(file, *dset_name):
	'''
	Read h5 file

	--- INPUT ---
	file        input h5 file
	dset_name   labels of data to read
	--- OUTPUT ---
	dset_data   data
	'''
	hf = H5.File(file+h5ext, 'r')
	dset_data = []
	for name in dset_name:
		data = hf.get(name)
		flag = 1
		if flag==1:
			data = np.array(data)
		dset_data.append(data)

	hf.close()
	data = np.array(data)

	return dset_data

def write_hdf5(data, name, file, append=False):
	'''
	Write h5 file

	--- INPUT ---
	data        data
	name        label of data
	file        file name of the new h5 file
	append      True if not overwrite (default: False)
	--- OUTPUT ---
	new h5 file
	'''
	if append==True:
		hf = H5.File(file+h5ext, 'a')
	else:
		hf = H5.File(file+h5ext, 'w')
	
	hf.create_dataset(name, data=data)

	hf.flush()
	hf.close()

def read_ascii(file):
	'''
	Read ASCII file

	--- INPUT ---
	file        input ASCII file
	--- OUTPUT ---
	'''
	f = open(file+ascext, 'r')
	## f.read() -> str | f.readlines() -> list
	data =[]
	for line in f.readlines():
		line = line.strip()
		print(line)
		if line[0]!='#':
			line = line.split()
			element = []
			for vec in line:
				element.append(vec)
			data.append(element)

	f.close()
	data = np.array(data)

	return data

def read_csv(file):
	'''
	Read csv file

	--- INPUT ---
	file        input csv file
	--- OUTPUT ---
	data        data
	'''
	with open(file+csvext, 'r') as csvfile:
		reader = csv.reader(csvfile)
		data = []
		for line in reader:
			data.append(line)

	return np.array(data)

def write_csv(file, header, data):
	'''
	Read fits file

	--- INPUT ---
	file        file name of the new csv file
	header      label of data
	data        data
	--- OUTPUT ---
	new csv file
	'''
	with open(file+csvext, 'w') as csvfile:
		writer = csv.writer(csvfile)
		writer.writerow(header)
		writer.writerows(data)
