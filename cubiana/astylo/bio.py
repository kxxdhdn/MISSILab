#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""

Basic Input & Output

"""

import numpy as np
from astropy.io import fits
import h5py as H5
import csv

global fitsext, h5ext
fitsext = '.fits'
h5ext = '.h5'
ascext = '.txt'
csvext = '.csv'

def read_fits(file, file_unc=None, wmod=0):
	'''
	Read fits file (auto detect dim)

	------ INPUT ------
	file                input fits file
	file_unc            input uncertainty file
	wmod                output wave mode (Default: 0)
	                      0 - 1darray; 
	                      1 - FITS_rec.
	------ OUTPUT ------
	cl                  output object
	  header              header of primary HDU
	  data                data in primary HDU
	  header_w            header of W-TAB
	  wave                data in table 1 (None if does not exist)
	  unc                 uncertainty array
	'''
	## Initialize output object
	cl = type('', (), {})()
	cl.header_w = None
	cl.wave = None
	cl.unc = None

	## Read header & data
	with fits.open(file+fitsext) as hdul:
		cl.HDUL = hdul
		hdr = hdul[0].header
		cl.data = hdul[0].data
		cl.header = hdr

		## Read wavelength
		if len(hdul)==2:
			cl.header_w = hdul[1].header
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
			
			cl.wave = wave
	
	if file_unc is not None:
		## Read uncertainty data
		with fits.open(file+fitsext) as hdul:
			cl.unc = hdul[0].data

	return cl

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

def read_hdf5(file, *header):
	'''
	Read h5 file

	------ INPUT ------
	file                input h5 file
	header              labels of data to read
	------ OUTPUT ------
	dataset             data
	'''
	hf = H5.File(file+h5ext, 'r')
	dataset = []
	for hdr in header:
		data = hf.get(hdr)
		data = np.array(data)
		dataset.append(data)

	hf.close()

	return dataset

def write_hdf5(file, header, data, append=False):
	'''
	Write h5 file

	------ INPUT ------
	file                file name of the new h5 file
	header              label of data (one at a time)
	data                data
	append              True: if not overwrite (default: False)
	------ OUTPUT ------
	'''
	if append==True:
		hf = H5.File(file+h5ext, 'a')
	else:
		hf = H5.File(file+h5ext, 'w')
	
	hf.create_dataset(header, data=data)

	hf.flush()
	hf.close()

def read_ascii(file, dtype=str, ascext=ascext):
	'''
	Read ASCII file

	------ INPUT ------
	file                input ASCII file
	dtype               data type (default: 'str')
	------ OUTPUT ------
	dataset             output data array
	'''
	with open(file+ascext, 'r') as f:
		## f.read() -> str | f.readlines() -> list
		dataset =[]
		for line in f.readlines():
			line = line.strip()
			# print(line)
			if line[0]!='#':
				line = list(map(dtype, line.split()))
				data = []
				for vec in line:
					data.append(vec)
				dataset.append(data)

	dataset = np.array(dataset)

	return dataset

def read_csv(file, *header):
	'''
	Read csv file

	------ INPUT ------
	file                input csv file
	header              labels of data to read
	------ OUTPUT ------
	dataset             dataset
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

	return dataset

def write_csv(file, header, dataset, append=False):
	'''
	Read fits file

	------ INPUT ------
	file                file name of the csv file
	header              labels of data, list('label1', 'label2', ...)
	dataset             data, list([d11, d12, ...], [d21, d22, ...], ...)
	------ OUTPUT ------
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
