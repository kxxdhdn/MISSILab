#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""

SYNTHETIC PHOTOMETRY

"""

from os import system
from astropy.io import fits
import numpy as np
import h5py
import matplotlib.pyplot as plt
import subprocess as SP


def synthetic_photometry(w_spec, Fnu_spec0, filt_UTF8):

	## Write the input file
	filIN = './synthetic_photometry_input.h5'

	hf = h5py.File(filIN, 'w')
	filt_ASCII = [n.encode("ascii","ignore") for n in filt_UTF8]
	ds_filt = hf.create_dataset("Filter label",dtype="S"+str(len(filt_ASCII[0])), \
	data=filt_ASCII,shape=(len(filt_UTF8),))
	#print("---------filt_UTF8--------", type(filt_UTF8))
	#print("--------filt_ASCII--------", len(filt_ASCII[0]))
	ds_w = hf.create_dataset('Wavelength (microns)', data=w_spec)
	ds_Fnu = hf.create_dataset('Flux (x.Hz-1)', data=Fnu_spec0)
	ds_bool = hf.create_dataset('(docalib,dophot)', data=[1,1])
	hf.flush()
	hf.close()

	## Call the Fortran code
	SP.call('synthetic_photometry')
	#system('synthetic_photometry')

	## Read the output
	filOUT = './synthetic_photometry_output.h5'

	hf = h5py.File(filOUT, 'r')
	Fnu_filt = hf.get('Flux (x.Hz-1)')
	Fnu_filt = np.array(Fnu_filt)
	filt = hf.get('Filter label')
	wcen = hf.get('Central wavelength (microns)')
	wcen = np.array(wcen)
	vmat = hf.get('Covariance matrix')
	rmat = hf.get('Correlation matrix')
	smat = hf.get('Standard dviation matrix')
	hf.close()


	## Final cleaning
	#system('rm ' + filIN)
	#system('rm ' + filOUT)

	return wcen, Fnu_filt


"""
------------------------------ MAIN (test) ------------------------------
"""
if __name__ == "__main__":
	
	path = './n66_LL/'
	filename = 'n66_LL1_cube'
	## select fixed rectangle: lower left (xb, yb) & upper right (xt, yt)
	xb, xt = 32, 44
	yb, yt = 47, 60
	
	with fits.open(path+filename+'.fits') as hdul:
		data = hdul[0].data
		hdr = hdul[0].header
		NAXIS1 = hdr['NAXIS1']
		NAXIS2 = hdr['NAXIS2']
		NAXIS3 = hdr['NAXIS3']
		wvl = hdul[1].data[0][0][:,0]
	#	wvl = hdul[1].data # for rewitten header
	
	## fluxes in box 
	data1 = data[:, yb:yt, xb:xt].reshape(NAXIS3, ((xt-xb)*(yt-yb)))
	flux0 = np.nansum(data1, axis=1).reshape((NAXIS3,1,1))
	
	## test (ref)
	#data_r = data[:, yb:yt, xb:xt]
	#flux_r = np.nansum(data_r, axis=(1,2)).reshape((NAXIS3,1,1))
	#print(flux0-flux_r)
	## test (reshape error)
	#data_s = data[:, yb:yt, xb:xt].reshape(((xt-xb)*(yt-yb)), NAXIS3).swapaxes(0,1)
	#flux_s = np.nansum(data_s, axis=1).reshape((NAXIS3,1,1))
	#data_t = data[:, yb:yt, xb:xt].reshape(((xt-xb)*(yt-yb)), NAXIS3).transpose(1,0)
	#flux_t = np.nansum(data_t, axis=1).reshape((NAXIS3,1,1))
	#print(data1-flux_t)
	#print("------")
	#print(data_s.shape)#-flux_r)
	
	## photometry
	wcen, Fnu_filt = synthetic_photometry(wvl, flux0, ["MIPS1"])
	print(wcen, Fnu_filt)
