#!/usr/bin/env python
# -*- coding: utf-8 -*

"""
calib.py is a routine for the inter-calibration of IRS spectra. It provides 
the inter-calibration with IRAC (1,2,3,4), MIPS, WISE bands.
"""

from astropy.io import fits
from astropy import units as u
from astropy import wcs
import numpy as np
import matplotlib.pyplot as plt

fits_path = './data/'
data_filename = 'm83_sl1_SL1_cube'
unc_filename = 'm83_sl1_SL1_cube_unc'

x = 125
y = 128

hdul = fits.open(fits_path+data_filename+'.fits')
data = hdul[0].data
hdr = hdul[0].header
NAXIS3 = hdr['NAXIS3']
wvl = hdul[1].data[0][0][:,0]
hdulunc = fits.open(fits_path+unc_filename+'.fits')
unc = hdulunc[0].data

band = [7.5, 10] # [lam_min, lam_max]
dlam = wvl[1]-wvl[0]
fbin = []
for k in range(NAXIS3):
	if (wvl[k]>band[0] and wvl[k]<band[1]):
		fbin.append(data[k,y,x]*dlam)
fband = np.nansum(fbin)
print("Band flux integral = {} MJy.Hz/sr".format(fband))

"""
plot
"""
fig = plt.figure("PLOT")
ax = plt.subplot()
ax.errorbar(wvl, data[:,y,x], yerr=unc[:,y,x], c='k', ls='-', lw=.5, ecolor= 'r', label="data_filename")
ax.set_title("Spectrum")
ax.set_xlabel("WVL (micron)")
ax.set_ylabel(r"$F_{\nu} (MJy/sr)$")
#ax.set_xscale('log', nonposy='clip')
#ax.set_yscale('log', nonposy='clip')
ax.legend(loc='upper left')

plt.show()
