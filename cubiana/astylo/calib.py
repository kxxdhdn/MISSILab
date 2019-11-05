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
from sinout import read_hdf5, write_hdf5, read_fits, read_ascii
from splot import plot2d_m, mycolorlib
from processim import iconvolve, project


##-----------------------------------------------

##			"PSFcalib" based tools

##-----------------------------------------------

class PSFcalib:
	'''
	PSF calibration (See Kospal08)
	'''
	def __init__(self):

		pass

##-----------------------------------------------

##			"intercalib" based tools

##-----------------------------------------------

class intercalib:
	'''
	Intercalibration
	'''
	def __init__(self, filIN=None, wmod=0):
		
		## INPUTS
		self.filIN = filIN

		if filIN is not None:
			self.hdr, self.im, self.wvl = read_fits(filIN, wmod=wmod)

	def synthetic_photometry(self, filt_UTF8, \
		w_spec=None, Fnu_spec=None, verbose=False):
		'''
		External Fortran library needed

		--- INPUT ---
		filt_UTF8    name list of photometry
		w_spec       wavelengths of obs. spectrum (via filIN)
		Fnu_spec     fluxes of spectrum (via filIN)
		verbose      keep i/o h5 file (Default: False)
		--- OUTPUT ---
		wcen         center wavelength
		Fnu_filt     fluxes of synthetic photometry
		std_dev      standard deviation
		'''
		if self.filIN is not None:
			n1 = self.hdr['NAXIS1']
			n2 = self.hdr['NAXIS2']
			n3 = self.hdr['NAXIS3']
			## Interpolation correction (spectrum debut)
			##-------------------------------------------
			w_spec = [.1]
			w_spec.extend(self.wvl) # extended wave
			Fnu_spec = np.zeros((n3+1)*n2*n1).reshape(n3+1, n2, n1)
			Fnu_spec[1:,:,:] = np.copy(self.im) # extended cube

		## Write the input file
		fortIN = './synthetic_photometry_input'

		filt_ASCII = [n.encode('ascii', 'ignore') for n in filt_UTF8]
		write_hdf5(fortIN, 'Filter label', filt_ASCII)
		write_hdf5(fortIN, 'Wavelength (microns)', w_spec, append=True)
		write_hdf5(fortIN, 'Flux (x.Hz-1)', Fnu_spec, append=True)
		write_hdf5(fortIN, '(docalib,dophot)', [1,1], append=True)

		## Call the Fortran lib
		##----------------------
		SP.call('synthetic_photometry')

		## Read the output
		fortOUT = './synthetic_photometry_output'

		wcen, Fnu_filt, smat = read_hdf5(fortOUT, \
			'Central wavelength (microns)', \
			'Flux (x.Hz-1)', \
			'Standard deviation matrix')
		wcen = wcen[0]
		Fnu_filt = Fnu_filt[0]
		std_dev = smat[0][0]

		## Clean temperary file
		##----------------------
		if verbose==False:
			SP.call(['rm', '-rf', fortIN+'.h5'])
			SP.call(['rm', '-rf', fortOUT+'.h5'])

		return wcen, Fnu_filt, std_dev
		
class spec2phot(intercalib):
	'''
	Intercalibration between spectrometry and photometry (REF)
	'''
	def __init__(self, filIN, filREF, wmod=0, ):
		super().__init__(filIN, wmod)

		pass

class phot2phot:
	'''
	Intercalibration between two photometry
	'''
	def __init__(self, filIN, filREF, filKER=None, saveKER=None, \
		wmod=0, uncIN=None, wmod_unc=0, Nmc=0, filOUT=None):

		## Convolution (optional)
		if filKER is not None:
			conv = iconvolve(filIN, filKER, saveKER, \
				wmod, uncIN, wmod_unc, filOUT=filPRO)
		else:
			filPRO = filIN

		## Reprojection
		self.pro = project(filPRO, filREF, filOUT=filOUT)

	def image(self):
		return self.pro.image()

def specorrect(filIN, factor=1., offset=0., wlim=(None,None), \
	wmod=0, filOUT=None):
	'''
	calibrate spectra from different obs. in order to eliminate gaps
	

	--- INPUT ---
	filIN       input fits file
	factor      scalar or ndarray (Default: 1.)
	offset      scalar or ndarray (Default: 0.)
	wlim        wave limits (Default: (None,None) = all spectrum)
	wmod        wave mode
	filOUT      overwrite fits file (Default: NO)
	--- OUTPUT ---
	im          im = factor * im + offset
	'''
	hdr, im, wvl = read_fits(file=filIN, wmod=wmod)

	if wlim[0] is None:
		wmin = wvl[0]
	else:
		wmin = wlim[0]
	if wlim[1] is None:
		wmax = wvl[-1]
	else:
		wmin = wlim[1]
		
	for k, lam in enumerate(wvl):
		if lam>=wmin and lam<=wmax:
			im[k,:,:] = factor * im[k,:,:] + offset
				
	if filOUT is not None:
		write_fits(filOUT, hdr, im, wave=wvl)

	return im


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
	greylines = []
	pinklines = []
	greylines.extend([2.3567863, 5.1532226]) # AKARI/IRC
	greylines.extend([5.242817, 7.597705]) # Spitzer/IRS-SL2
	pinklines.extend([7.3675313, 8.66892]) # Spitzer/IRS-SL3
	greylines.extend([7.5337057, 14.736635]) # Spitzer/IRS-SL1
	greylines.extend([14.266611, 21.051888]) # Spitzer/IRS-LL2
	pinklines.extend([19.483675, 21.50092]) # Spitzer/IRS-LL3
	greylines.extend([20.555237, 38.41488]) # Spitzer/IRS-LL1
	p.ax.vlines(greylines, 0, 1.1, linestyles='dotted', colors='grey')
	p.ax.vlines(pinklines, 0, 1.1, linestyles='dotted', colors='pink')

	## tick setting
	##-------------------- x --------------------------
	xtic = [2, 3, 4, 5, 6, 7, 8, 10, 12, 15, 20, 30, 40]
	xtic_min = np.arange(2., 41., 1.)
	p.ax.set_xticks(xtic, minor=False) # major
	p.ax.set_xticks(xtic_min, minor=True) # minor
	# ScalarFormatter().set_scientific(False)
	p.ax.xaxis.set_major_formatter(ScalarFormatter()) # major
	p.ax.xaxis.set_minor_formatter(NullFormatter()) # minor
	# p.ax.minorticks_off()
	##--------------------- y --------------------------
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
	
	##--------------------
	## photometry_profile
	##--------------------
	# p = photometry_profile('./data/', \
	# 	'IRAC1', 'IRAC2', 'IRAC3', 'IRAC4', 'MIPS1', \
	# 	'WISE1', 'WISE2', 'WISE3', 'WISE4')
	# p.show()

	##---------------------------------
	## intercalib.synthetic_photometry
	##---------------------------------
	calib = intercalib('/Users/dhu/data/pahpedia/M83/output/M83')
	wcen, Fnu, sig = calib.synthetic_photometry(['MIPS1'])
