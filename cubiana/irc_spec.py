#!/usr/bin/env python
# -*- coding: utf-8 -*-

import logging, sys
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
# from astylo.calib import intercalib, phot2phot
# from astylo.plot import plot2d

## Local
from param import (
	src, Nmc, path_cur, root, path_irc, fits_irc, parobs, 
	path_out, path_tmp, path_build, path_tests, verbose, 
)
##---------------------------
##       Initialisation
##---------------------------

##---------------------------
##        Build slits
##---------------------------
print(parobs)
