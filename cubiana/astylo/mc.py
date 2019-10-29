#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""

Uncertainty propagation with Monte-Carlo method

"""

import numpy as np

## astylo
from sinout import write_fits


class MCerror:
	'''
	Monte-Carlo propagated error calculator
	'''
	pass

def MunC(arr, axis=0, filOUT=None, hdr=None, wave=None):
	'''
	Monte-Carlo propagated uncertainties 
	
	--- INPUT ---
	arr         arr
	axis        axis or axes along which unc is calculated
	filOUT      Default: no output fits file
	hdr         fits file header
	wave        wavelengths if 3D
	--- OUTPUT ---
	err         one dim less than arr
	'''
	## Detect dimension
	##------------------
	ax = arr.shape
	NAXIS = np.size(ax)
	dim = NAXIS - 1
	if axis==0:
		if dim==2:
			err = np.full((ax[1],ax[2]), np.nan)
			for i in range(ax[2]):
				for j in range(ax[1]):
					err[j,i] = np.nanstd(arr[:,j:i], ddof=1)
		elif dim==3:
			err = np.full((ax[1],ax[2],ax[3]), np.nan)
			for i in range(ax[3]):
				for j in range(ax[2]):
					for k in range(ax[1]):
						err[k,j,i] = np.nanstd(arr[:,k,j,i], ddof=1)
		else:
			print('ERROR: dimension not supported! ')
	else:
		print('ERROR: array format not supported! ')

	if filOUT!=None:
		comment = "Monte-Carlo propagated uncertainty file."
		write_fits(filOUT, header=hdr, data=err, wave=wave, COMMENT=comment)
	
	return err
