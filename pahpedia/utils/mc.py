#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
from utils.rwfits import *

def calunc(data, NAXIS, filOUT=None, hdr=None, wvl=None):
	
	## detect dimension
	dim = np.size(NAXIS)
	if dim==1:
		err = np.full(NAXIS[0], np.nan)
		for i in range(NAXIS[0]):
			err[i] = np.nanstd(data[:,i], ddof=1)
	elif dim==2:
		err = np.full((NAXIS[1], NAXIS[0]), np.nan)
		for i in range(NAXIS[0]):
			for j in range(NAXIS[1]):
				err[j,i] = np.nanstd(data[:,j,i], ddof=1)
	else:
		err = np.full((NAXIS[2], NAXIS[1], NAXIS[0]), np.nan)
		for i in range(NAXIS[0]):
			for j in range(NAXIS[1]):
				for k in range(NAXIS[2]):
					err[k,j,i] = np.nanstd(data[:,k,j,i], ddof=1)

	if filOUT!=None:
		comment = "Monte-Carlo propagated uncertainty file."
		write_fits(filOUT, data=err, hdr=hdr, wvl=wvl, COMMENT=comment)
	
	return err
