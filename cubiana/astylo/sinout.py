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

def write_fits(file, header, data, wave=None, **hdrl):
	'''
	Write fits file

	--- INPUT ---
	file        input fits file
	header      header of primary HDU
	data        data in primary HDU
	wave        data in table 1 (default: None)
	--- OUTPUT ---
	new fits file
	'''
	for key, value in hdrl.items():
		header[key] = value
	primary_hdu = fits.PrimaryHDU(header=header, data=data)
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
	header      header of primary HDU
	w           2D WCS
	is3d        if input data is 3D: True
	'''
	hdr = fits.open(file+fitsext)[0].header

	header = hdr.copy()
	if header['NAXIS']==3:
		is3d = True
		for kw in hdr.keys():
			if '3' in kw:
				del header[kw]
		header['NAXIS'] = 2
		header['COMMENT'] = "This header is adapted to 2D WCS extraction need. "
	else:
		is3d = False
	
	w = WCS(header, naxis=2)

	return header, w, is3d

def read_hdf5(file, *header):
	'''
	Read h5 file

	--- INPUT ---
	file        input h5 file
	header      labels of data to read
	--- OUTPUT ---
	dataset     data
	'''
	hf = H5.File(file+h5ext, 'r')
	dataset = []
	for hdr in header:
		data = hf.get(hdr)
		data = np.array(data)
		dataset.append(data)

	hf.close()

	return np.array(dataset)

def write_hdf5(file, header, data, append=False):
	'''
	Write h5 file

	--- INPUT ---
	file        file name of the new h5 file
	header      label of data (one at a time)
	data        data
	append      True if not overwrite (default: False)
	--- OUTPUT ---
	new h5 file
	'''
	if append==True:
		hf = H5.File(file+h5ext, 'a')
	else:
		hf = H5.File(file+h5ext, 'w')
	
	hf.create_dataset(header, data=data)

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
	dataset =[]
	for line in f.readlines():
		line = line.strip()
		print(line)
		if line[0]!='#':
			line = line.split()
			data = []
			for vec in line:
				data.append(vec)
			dataset.append(data)

	f.close()
	dataset = np.array(dataset)

	return dataset

def read_csv(file, *header):
	'''
	Read csv file

	--- INPUT ---
	file        input csv file
	header      labels of data to read
	--- OUTPUT ---
	dataset     dataset
	'''
	with open(file+csvext, 'r', newline='') as csvfile:
		reader = csv.DictReader(csvfile)
		dataset = []
		for hdr in header:
			data = []
			for row in reader:
				data.append(row[hdr])
			data = np.array(data)
			dataset.append(data)

	return np.array(dataset)

def write_csv(file, header, dataset, append=False):
	'''
	Read fits file

	--- INPUT ---
	file        file name of the csv file
	header      labels of data, list('label1', 'label2', ...)
	dataset     data, list([d11, d12, ...], [d21, d22, ...], ...)
	--- OUTPUT ---
	new csv file
	'''
	if append==True:
		mod = 'a'
	else:
		mod = 'w'

	with open(file+csvext, mod, newline='') as csvfile:
		writer = csv.DictWriter(csvfile, fieldnames=header)

		writer.writeheader()

		for i in range(len(dataset)):
			## Init dict
			row = {hdr: [] for hdr in header}
			data = dataset[i]
			## Write dict
			for j in range(len(header)):
				row[header[j]] = data[j]
			## Write csv row
			writer.writerow(row)
