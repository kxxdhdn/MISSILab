#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""

This is the visualization of the test(_profiles) of auxil.f90

"""
import os
import numpy as np
import matplotlib.pyplot as plt

## astylo
from astylo.iolib import read_hdf5
from astylo.plib import plot2d_m, plot2d

## Path
mroot = os.path.dirname(os.path.abspath(__file__))+'/'
path_out = mroot+'out/'
filename = 'test_profiles'

## Read h5 file
wvl = read_hdf5(path_out+filename, 'Wavelength (microns)')
nu = read_hdf5(path_out+filename, 'Wavelength (Hz)')
FnuBB = read_hdf5(path_out+filename, 'FnuBB (MJyovsr)')
FnuLINE = read_hdf5(path_out+filename, 'FnuLINE (MJyovsr)')
FnuBAND = read_hdf5(path_out+filename, 'FnuBAND (MJyovsr)')
Fnu = read_hdf5(path_out+filename, 'Fnu (MJyovsr)')
y = read_hdf5(path_out+filename, 'FnuTEST (MJyovsr)')
x = read_hdf5(path_out+filename, 'x')

## Normalization test
INTest = np.sum((x[1]-x[0])*y) # Trapezoidal rule
print('I(test) = ', INTest)

# plt.plot(wvl, FnuBB)
# plt.xlim(2., 40.)
# plt.show()
# exit()

p = plot2d_m(wvl, [Fnu, FnuBB, FnuLINE, FnuBAND], xall=wvl, 
	         xlog=1, ylog=1, xlim=(2., 40.), 
             xlab=r'$\lambda\,(micron)$', ylab=r'$F_{\nu} \,(MJy/sr)$', 
             lablist=['Total', 'Black Body', 'Gaussian', 'Split Lorentzian'],
             legend='best', cl=['k', 'orange', 'g', 'b'])
p.save(path_out+filename+'.png')
p.show()

print('>>> Coucou test_profiles [Done] <<<')
