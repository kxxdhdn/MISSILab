#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""

SYNTHETIC PHOTOMETRY

"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter, NullFormatter
import subprocess as SP

## astylo
from sinout import read_hdf5, write_hdf5, read_ascii
from splot import plot2d_m, mycolorlib

def synthetic_photometry(w_spec, Fnu_spec0, filt_UTF8):
	'''
	Extra Fortran library needed

	--- INPUT ---
	w_spec       wavelengths of obs. spectrum
	Fnu_spec0    fluxes of spectrum
	filt_UTF8    name list of photometry
	--- OUTPUT ---
	wcen         center wavelength
	Fnu_filt     fluxes of synthetic photometry
	std_dev      standard deviation
	'''
	## Write the input file
	filIN = './synthetic_photometry_input'

	filt_ASCII = [n.encode('ascii', 'ignore') for n in filt_UTF8]
	write_hdf5(filIN, 'Filter label', filt_ASCII)
	write_hdf5(filIN, 'Wavelength (microns)', w_spec, append=True)
	write_hdf5(filIN, 'Flux (x.Hz-1)', Fnu_spec0, append=True)
	write_hdf5(filIN, '(docalib,dophot)', [1,1], append=True)

	## Call the Fortran lib
	##----------------------
	SP.call('synthetic_photometry')

	## Read the output
	filOUT = './synthetic_photometry_output'

	wcen, Fnu_filt, smat = read_hdf5(filOUT, \
		'Central wavelength (microns)', \
		'Flux (x.Hz-1)', \
		'Standard deviation matrix')
	wcen = wcen[0]
	Fnu_filt = Fnu_filt[0]
	std_dev = smat[0][0]

	## Clean temperary file
	##----------------------
	SP.call(['rm', '-rf', filIN+'.h5'])
	SP.call(['rm', '-rf', filOUT+'.h5'])

	return wcen, Fnu_filt, std_dev

def photometry_profile(path_dat, *photometry):
	'''
	--- INPUT ---
	path_dat     profile data path
	fig_save     save profiles
	photometry   photometry
	--- OUTPUT ---
	'''
	## Read data
	##-----------
	lam = []
	val = []
	for phot in photometry:
		dat = read_ascii(path_dat+phot, dtype=float)
		lam.append(dat[:,0])
		val.append(dat[:,1])
	lam = np.array(lam)
	val = np.array(val)

	## Plotting setting
	##------------------
	p = plot2d_m(lam, val, xlim=(1.9, 40.), ylim=(-.01, 1.01), \
	xlog=1, ylog=1, cl=mycolorlib[2:], lw=1.8, \
	lablist=photometry, \
	xlab=r'$Wavelength,\,\,\lambda\,\,[\mu m]$', \
	ylab='Response', \
	# ylab='Spectral response\n(electrons / photon)', \
	legend='upper left', figsize=(12,3))

	p.set_border(left=.05, bottom=.2, right=.99, top=.99)

	# sizeXL = 50
	# p.set_font(xticksize=sizeXL, yticksize=sizeXL, \
	# 	axesize=sizeXL, legendsize=sizeXL)

	## vlines (e.g. band markers)
	##----------------------------
	'''
	AKARI/IRC: 2.3567863, 5.1532226
	Spitzer/IRS-SL2: 5.242817, 7.597705
	Spitzer/IRS-SL1: 7.5337057, 14.736635
	Spitzer/IRS-LL2: 14.266611, 21.051888
	Spitzer/IRS-LL1: 20.555237, 38.41488
	'''
	lines = [2.3567863, 5.1532226, \
	5.242817, 7.597705, \
	7.5337057, 14.736635, \
	14.266611, 21.051888, \
	20.555237, 38.41488]
	p.ax.vlines(lines, 0, 1.1, linestyles='dotted', colors='grey')

	## tick setting
	##-----------------------------------------------
	xtic = [2, 3, 4, 5, 6, 7, 8, 10, 12, 15, 20, 30, 40]
	xtic_min = np.arange(2., 41., 1.)
	p.ax.set_xticks(xtic, minor=False) # major
	p.ax.set_xticks(xtic_min, minor=True) # minor
	# ScalarFormatter().set_scientific(False)
	p.ax.xaxis.set_major_formatter(ScalarFormatter()) # major
	p.ax.xaxis.set_minor_formatter(NullFormatter()) # minor
	# p.ax.minorticks_off()
	##-----------------------------------------------
	ytic = np.arange(0, 1.01, .2)
	ytic_min = np.arange(0, 1., .1)
	p.ax.set_yticks(ytic, minor=False) # major
	p.ax.set_yticks(ytic_min, minor=True) # minor
	# ScalarFormatter().set_scientific(False)
	p.ax.yaxis.set_major_formatter(ScalarFormatter()) # major
	p.ax.yaxis.set_minor_formatter(NullFormatter()) # minor
	# p.ax.minorticks_off()
	##-----------------------------------------------

	return p

"""
------------------------------ MAIN (test) ------------------------------
"""
if __name__ == "__main__":
	
	p = photometry_profile('./data/', \
		'IRAC1', 'IRAC2', 'IRAC3', 'IRAC4', 'MIPS1', \
		'WISE1', 'WISE2', 'WISE3', 'WISE4')

	# file_save = ''
	# p.save(file_save, transparent=True)
	p.show()
