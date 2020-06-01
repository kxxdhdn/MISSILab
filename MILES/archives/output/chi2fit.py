#!/usr/bin/env python
# -*- coding: utf-8 -*-

from sinout import read_hdf5
from splot import plot2d, plot2d_multi
import numpy as np
import matplotlib.pyplot as plt

## read_optics
path = '/Users/dhu/Dropbox/MILES/output/'
filename = path+'chi2fit'
wvl, Fnu, FnuBB, FnuLINE, FnuBAND, FnuSTAR = read_hdf5(filename, \
	'Wavelength (micron)', 'Fnu (MJyovsr)', \
	'FnuBB (MJyovsr)', 'FnuLINE (MJyovsr)', 'FnuBAND (MJyovsr)', \
	'FnuSTAR (MJyovsr)')

plot2d_multi(wvl, [Fnu, FnuBB, FnuLINE, FnuBAND, FnuSTAR], xlog=1, ylog=1, \
	xlab=r'$\lambda\,(micron)$' , \
	ylab=r'$F_{\nu} \,(MJy/sr)$', \
	lab=['Total', 'Black Body', 'Gaussian', 'Split Lorentzian', 'Star'], \
	savenm=path+'test_fit.png')

plt.show()
