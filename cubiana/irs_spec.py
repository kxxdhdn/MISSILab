#!/usr/bin/env python
# -*- coding: utf-8 -*-

import logging, sys
# logging.basicConfig(stream=sys.stderr, level=logging.ERROR)
# print(logging.getLogger())
logging.disable(sys.maxsize)
import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning) 

from tqdm import tqdm, trange

import os
import numpy as np

## astylo
from astylo.bio import read_fits, write_fits, read_ascii
from astylo.proc import islice, iconvolve, iswarp, wclean
from astylo.lib import fclean, pix2sr
from astylo.plot import plot2d

## Local
from param import (
	src, Nmc, path_cur, path_idl, 
	root, path_irs, path_phot, path_ker, 
	fits_irs, chnl, fits_ker, path_out, csv_ker, 
	phot, phot0, path_cal, fits_phot, fits_phot0, 
	path_tmp, path_slices, path_tests, verbose, 
)

##---------------------------
##      Initialisation
##---------------------------
Nch = 4 # Number of chnl
'''
##---------------------------
##       Combine obs
##---------------------------
swp = iswarp(sum(fits_irs, []), \
	center='9:55:52,69:40:45', pixscale='1.67', \
	verbose=False, tmpdir=path_tmp)

for i in trange(Nch, #leave=False, 
	desc='<iswarp> IRS Combining ({} chnl)'.format(Nch)):
	for j in trange(Nmc+1, leave=False, 
		desc='IRS Combining [MC]'):
		if j==0:
			comb = swp.combine(fits_irs[i], \
				'wgt_avg', keepedge=True, \
				tmpdir=path_tmp+'MC_no/', \
				filOUT=path_tmp+src+'_'+chnl[i])
		else:
			comb = swp.combine(fits_irs[i], 'wgt_avg', \
				keepedge=True, uncpdf='norm', \
				tmpdir=path_tmp+'MC_'+str(j)+'/', \
				filOUT=path_tmp+src+'_'+chnl[i]+'_'+str(j))

## PSF Convolution
##-----------------
for i in trange(Nch, #leave=False, 
	desc='<iconvolve> IRS Smoothing ({} chnl)'.format(Nch)):
	for j in trange(Nmc+1, leave=False, 
		desc='<iconvolve> IRS Smoothing [MC]'):
		if j==0:
			conv = iconvolve(path_tmp+src+'_'+chnl[i], \
				fits_ker, csv_ker, \
				filTMP=path_slices+src+'_'+chnl[i], \
				filOUT=path_out+src+'_'+chnl[i])
		else:
			conv = iconvolve(path_tmp+src+'_'+chnl[i]+'_'+str(j), \
				fits_ker, csv_ker, \
				filTMP=path_slices+src+'_'+chnl[i]+'_'+str(j), \
				filOUT=path_tmp+src+'_'+chnl[i]+'_conv_'+str(j))
		conv.do_conv(ipath=path_idl)
'''
for i in trange(Nch, #leave=False, 
	desc='IRS Cal unc ({} chnl)'.format(Nch)):
	mcimage = []
	for j in trange(Nmc+1, leave=False, 
		desc='IRS Reading [MC]'):
		if j==0:
			hd0 = read_fits(path_out+src+'_'+chnl[i])
			header = hd0.header
			wvl = hd0.wave
		else:
			hd = read_fits(path_tmp+src+'_'+chnl[i]+'_conv_'+str(j))
			mcimage.append(hd.data)
	mcimage = np.array(mcimage)
	unc = np.nanstd(mcimage, axis=0)
	write_fits(path_out+src+'_'+chnl[i]+'_unc', header, unc, wvl)
