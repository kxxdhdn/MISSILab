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
filename = 'test_gsamp'

## Read h5 files
wvl, FnuOBS, FnuMOD = \
    read_hdf5(path_out+filename, 
    'Wavelength (microns)', 'FnuOBS (MJyovsr)', 'FnuMOD (MJyovsr)')
Fnu1par = read_hdf5(path_out+filename, 
    'var1par')
# print(Fnu1par[0].shape)
# exit()

p = plot2d_m(wvl, Fnu1par[0], 
             xall=wvl, xlog=1, ylog=1, xlim=(1., 40.), ylim=(0,2.), 
             xlab=r'$\lambda\,(micron)$', ylab=r'$F_{\nu} \,(MJy/sr)$', 
             lablist=[str(i) for i in range(5)], 
             legend='best', cl=['r', 'orange', 'y', 'g', 'b', 'purple'])
# p = plot2d_m(wvl, [FnuOBS, FnuMOD], 
#              xall=wvl, xlog=1, ylog=1, xlim=(2., 40.), ylim=(0,2.), 
#              xlab=r'$\lambda\,(micron)$', ylab=r'$F_{\nu} \,(MJy/sr)$', 
#              lablist=['Obs', 'Mod'], 
#              legend='best', cl=['y', 'k', 'pink', 'g', 'b', 'orange'])
p.save(path_out+filename+'.png')
p.show()

print('>>> Coucou test_gsamp [Done] <<<')
