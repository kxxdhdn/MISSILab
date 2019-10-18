#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""

This is a simple test of read_optics.f90 and 
blackbody() function in statistics.f90

"""

from sinout import read_ascii, read_hdf5
from splot import plot2d, plot2d_multi
import numpy as np
import matplotlib.pyplot as plt
from astropy.modeling.blackbody import blackbody_nu

## read_optics
path = '/Users/dhu/Dropbox/MILES/output/'
filename = path+'test_modifBB'
data = read_hdf5(filename, 'Wavelength (micron)', \
	'Wavelength (Hz)', 'Fnu (MJyovsr)')
# data = read_ascii('test_modifBB.txt', Nvec=2)[0]
wvl = data[0]
nu = data[1]

## astropy
fopt = '/Users/dhu/PhD/SwING/Model_templates/Cross_sections/Data/'\
	+'optics_Sil_D03'
radius, rho, Qabs = read_hdf5(fopt, \
	'Grain radius (microns)', 'Grain density (kg.m-3)', \
	'Grain absorption efficiency (Qabs)')
a = radius[351]
Qabs = np.array(Qabs)[:,351]
dist = 1.
temp = 100.
mass = 3.
Msun = 1.989E30
kpc = 3.086E19
Jy = 1.E-26
Fnu = mass/dist**2 * Msun/kpc**2 \
	* 3.*np.pi/4./rho * Qabs/a \
	* blackbody_nu(nu, temp).value \
	* 1.E-7 * 1.E4 / 1.E6 # erg/s/cm−2/sr/Hz -> MW/m2/sr/Hz = MJy/sr

## plot
plot2d_multi(wvl, [data[2], Fnu], xlog=1, ylog=1, \
	xlab=r'$\lambda\,(micron)$' , \
	ylab=r'$F_{\nu}/\nu \,(MJy/sr/Hz)$', \
	lab=['read_optics','astropy'], \
	savenm=path+'Sil_3Msun_100K.png')

plt.show()
