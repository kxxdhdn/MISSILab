#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""

This is the visualization of the test(_profiles) of auxil.f90

"""
import os
import numpy as np
import matplotlib.pyplot as plt

## rapyuta
from rapyuta.inout import read_hdf5
from rapyuta.plots import pplot

## Path
mroot = os.path.dirname(os.path.abspath(__file__))+'/'
path_out = mroot+'out/'
filename = 'test_profiles'

## Read h5 file
wvl = read_hdf5(path_out+filename, 'Wavelength (microns)')
nu = read_hdf5(path_out+filename, 'Wavelength (Hz)')
FnuCONT = read_hdf5(path_out+filename, 'FnuCONT (MJyovsr)')
FnuLINE = read_hdf5(path_out+filename, 'FnuLINE (MJyovsr)')
FnuBAND = read_hdf5(path_out+filename, 'FnuBAND (MJyovsr)')
Fnu = read_hdf5(path_out+filename, 'Fnu (MJyovsr)')
y = read_hdf5(path_out+filename, 'FnuTEST (MJyovsr)')
x = read_hdf5(path_out+filename, 'x')

## Normalization test
INTest = np.sum((x[1]-x[0])*y) # Trapezoidal rule
print('I(test) = ', INTest)

p = pplot(wvl, Fnu,
          xlog=1, ylog=1, xlim=(2., 40.), ylim=(1e-4,1e2),
          xlabel=r'$\lambda\,(micron)$',
          ylabel=r'$F_{\nu} \,(MJy/sr)$', 
          label='Total', loc='upper left',
          clib=['k', 'orange', 'g', 'b'])
p.add_plot(wvl, FnuCONT, label='Black Body')
p.add_plot(wvl, FnuLINE, label='Gaussian')
p.add_plot(wvl, FnuBAND, label='Split Lorentzian')

p.save(path_out+filename+'.png')
# p.show()

print('>>> Coucou test_profiles [Done] <<<')
