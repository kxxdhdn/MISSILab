#!/usr/bin/env python
# -*- coding: utf-8 -*-

import logging, sys
# logging.basicConfig(stream=sys.stderr, level=logging.ERROR)
# print(logging.getLogger())
logging.disable(sys.maxsize)

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
	fits_irs, fits_unc_irs, chnl, fits_ker, path_out, 
	phot, phot0, path_cal, fits_phot, fits_phot0, 
	path_test, path_tmp, verbose, 
)

##---------------------------
##      Initialisation
##---------------------------
## Set output path
fits_out_irs = path_out+src+'_IRS'
fits_out_unc_irs = path_out+src+'_IRS_unc'

##---------------------------
##       Combine obs
##---------------------------
for j in trange(Nmc, desc='imontage : Combining IRS...'):
	## Reproject all chnl to the last one in the chnl list
	mont = imontage(fits_irs, fits_irs[0], fmod='ext', ext_pix=2)
	mont.make()
	mont.combine(fits_out_irs, 'wgt_avg', fits_unc_irs, Nmc=Nmc, write_mc=True)#, do_rep=False)

print('imontage : Combining IRS cube [done]')

##---------------------------
##         Path Store
##---------------------------
