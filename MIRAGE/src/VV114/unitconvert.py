#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""

DustPedia image unit conversion

"""

# import warnings
# warnings.filterwarnings("ignore", category=RuntimeWarning)


import numpy as np
from matplotlib.ticker import ScalarFormatter, NullFormatter

## rapyuta
from rapyuta.inout import fclean, fitsext, read_fits, write_fits, read_hdf5
from rapyuta.imaging import Jy_per_pix_to_MJy_per_sr

## Local
from buildinfo import ( src, path_cal, path_out, path_phot )


## Banner
print('\n============================================================\n')

print('        MIRAGE - DustPedia image unit conversion - '+src)

print('\n============================================================\n')


## Convert DustPedia photometry unit (Jy/pix -> MJy/sr)
##------------------------------------------------------
phot = ['IRAC4']
raw_phot = path_phot+src+'_'+phot[0]+'_DP'
out_phot = path_cal+'_'+src+'_IRAC4_DP'
Jy_per_pix_to_MJy_per_sr(raw_phot, filOUT=out_phot)
