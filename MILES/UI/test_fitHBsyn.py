#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""

This is the visualization of the test(_specModel) of auxil.f90

"""
import os
import math
import numpy as np
import matplotlib.pyplot as plt
# from astropy.constants import M_sun, pc, 
# from astropy.units import Jy, 

## astylo
from astylo.iolib import read_hdf5
from astylo.plib import plot2d_m, plot2d

## Path
mroot = os.path.dirname(os.path.abspath(__file__))+'/../'
path_out = mroot+'tests/out/'
filename = 'test_fitHBsyn'

## Read h5 files
wvl, FnuOBS, FnuMOD, FnuCONT, FnuLINE, FnuBAND, FnuSTAR = \
    read_hdf5(path_out+filename, 
    'Wavelength (microns)', 'FnuOBS (MJyovsr)', 'FnuMOD (MJyovsr)', 
    'FnuCONT (MJyovsr)', 'FnuLINE (MJyovsr)', 'FnuBAND (MJyovsr)', 
    'FnuSTAR (MJyovsr)')
# wvl, FnuOBS, meanpar, stdevpar = read_hdf5(path_out+filename, 
#     'Wavelength (microns)', 'FnuOBS (MJyovsr)', 
#     'Mean of parameter value', 'Sigma of parameter value')

# ## Extract par
# i0 = 0
# Ncont = 3
# Mcont = []
# Tcont = []
# for i in range(Ncont):
#     Mcont.append(meanpar[i0+2*i])
#     Tcont.append(meanpar[i0+2*i+1])

# i0 = i0 + 2*Ncont
# Nline = 12
# Iline = []
# Cline = []
# Wline = []
# for i in range(Nline):
#     Iline.append(meanpar[i0+3*i])
#     Cline.append(meanpar[i0+3*i+1])
#     Wline.append(meanpar[i0+3*i+2])

# i0 = i0 + 3*Nline
# Nband = 14
# Iband = []
# Cband = []
# WSband = []
# WLband = []
# for i in range(Nband):
#     Iband.append(meanpar[i0+4*i])
#     Iband.append(meanpar[i0+4*i+1])
#     WSband.append(meanpar[i0+4*i+2])
#     WLband.append(meanpar[i0+4*i+3])

# i0 = i0 + 4*Nband
# Av = meanpar[i0]

# i0 = i0 + 1
# Fstar = meanpar[i0]

# ## Calculate model
# for i in range(Ncont):
#     FnuCONT = Mcont[i] * M_sun/pc**2 * 
# FnuBAND = 0
# FnuSTAR = 0
# Pabs = 0
# FnuLINE = 0

# FnuMOD = (FnuCONT + FnuBAND + FnuSTAR) * Pabs + FnuLINE

p = plot2d_m(wvl, [FnuOBS, FnuMOD, FnuCONT, FnuLINE, FnuBAND, FnuSTAR], 
             xall=wvl, xlog=1, ylog=1, xlim=(2., 40.), ylim=(0,2.), 
             xlab=r'$\lambda\,(micron)$', ylab=r'$F_{\nu} \,(MJy/sr)$', 
             lablist=['Obs', 'Total', 'cont', 'lines', 'bands', 'star'], 
             legend='best', cl=['y', 'k', 'pink', 'g', 'b', 'orange'])
p.save(path_out+filename+'.png')
p.show()

print('>>> Coucou test_fitHBsyn [Done] <<<')
