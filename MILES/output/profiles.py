#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""

This is the visualization of the test(_profiles) of auxil.f90

"""

from sinout import read_hdf5
from splot import plot2d, plot2d_multi
import numpy as np
import matplotlib.pyplot as plt

## read_optics
path = '/Users/dhu/Dropbox/MILES/output/'
filename = path+'test_profiles'
wvl, nu, FnuBB, FnuLINE, FnuBAND, Fnu, y, x = read_hdf5(filename, \
	'Wavelength (micron)', 'Wavelength (Hz)', \
	'FnuBB (MJyovsr)', 'FnuLINE (MJyovsr)', 'FnuBAND (MJyovsr)', \
	'Fnu (MJyovsr)', 'FnuTEST (MJyovsr)', 'x')
## normalization test
INTest = np.sum((x[1]-x[0])*y) # Trapezoidal rule
print('Should be 1: ', INTest)

plot2d_multi(wvl, [Fnu, FnuBB, FnuLINE, FnuBAND], xlog=1, ylog=1, \
	xlab=r'$\lambda\,(micron)$' , \
	ylab=r'$F_{\nu} \,(MJy/sr)$', \
	lab=['Total', 'Black Body', 'Gaussian', 'Split Lorentzian'], \
	savenm=path+'test_profiles.png')

plt.show()
