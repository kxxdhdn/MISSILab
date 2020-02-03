#!/usr/bin/env python
# -*- coding: utf-8 -*-

import logging, sys
# logging.basicConfig(stream=sys.stderr, level=logging.ERROR)
# print(logging.getLogger())
logging.disable(sys.maxsize)
import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning) 

import time
t0 = time.time()
from tqdm import tqdm, trange

import os
import numpy as np

## astylo
from astylo.bio import read_fits, write_fits, read_ascii
from astylo.proc import icrop, iconvolve, imontage, wclean
from astylo.lib import fclean, pix2sr
from astylo.plot import plot2d

## Local
from param import (
	src, Nmc, path_cur, path_idl, 
	root, path_irs, path_phot, path_ker, 
	fits_irs, fits_irs_unc, chnl, 
	fits_ker, path_out, csv_ker, 
	phot, phot0, path_cal, fits_phot, fits_phot0, 
	path_test, path_tmp, verbose, 
)

##---------------------------
##      Initialisation
##---------------------------
Nch = 4 # Number of chnl used
## Set output path
fits_out_irs = []
for i in range(Nch):
	fits_out_irs.append(path_out+src+'_'+chnl[i])

##---------------------------
##       Combine obs
##---------------------------
## Reproject all chnl to the last one in the chnl list
##--------------------
## Make ref frame header
ref_irs = fits_irs[0][0] # [chnl][label]
# print(ref_irs)
# exit()
mont = imontage(sum(fits_irs,[]), ref_irs, \
	fmod='ext', ext_pix=2, ftmp=path_tmp)
mont.make()
mont.footprint(path_tmp+'footprint')
# from astropy import wcs
# print(wcs.WCS(mont.hdr_ref))
# print(wcs.WCS(read_fits(path_tmp+'footprint0').header))
# exit()
# mont.hdr_ref = read_fits(path_tmp+'footprint0').header
print(fits_irs)
exit()
combim = []
for i in trange(Nch, leave=False, \
	desc='Building IRS cubes'):
	mont.combine(fits_irs[i], fits_out_irs[i], 'wgt_avg', \
		fits_irs_unc[i], Nmc=Nmc, do_rep=True)

## PSF Convolution
##-----------------
# for j in trange(Nmc+1, desc='iconvolve: Convolving PSF'):
# 	for f in fits_irs:
# 		filename = os.path.basename(f)
# 		f_rep = mont.path_tmp + filename+'_rep/'
# 		if j==0:
# 			conv = iconvolve(f_rep+'rep', fits_ker, csv_ker)
# 		else:
# 			conv = iconvolve(f_rep+'rep_'+str(j))

# 		conv.do_conv(ipath=path_idl)

# print('iconvolve: Convolving PSF...[done]')

print('Building IRS cubes...[done]')
