#!/usr/bin/env python
# -*- coding: utf-8 -*-

# import logging, sys
# logging.disable(sys.maxsize)
# import warnings
# warnings.filterwarnings("ignore", category=RuntimeWarning) 

from tqdm import tqdm, trange

import os
import numpy as np
from scipy.optimize import curve_fit
from scipy.stats import gmean
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter, LogFormatter, NullFormatter
## astylo
from iolib import read_fits, write_fits
from proc import iconvolve, iswarp, sextract
from calib import intercalib
from mlib import fclean, f_lin, f_lin0
from alib import pix2sr, get_pc

## Local
from param import (
	src, Nmc, path_irc, fits_irc, parobs, 
	path_idl, path_ker, csv_ker, 
	phot, phot0, path_phot, path_cal, fits_phot, fits_phot0, 
	path_out, path_tmp, path_build, path_tests, verbose, 
)
##---------------------------
##       Initialisation
##---------------------------
Nobs = len(parobs)
Nmc = 2

##---------------------------
##        Build slits
##---------------------------
'''
for i in range(Nobs):
	sext = sextract(path_irc, parobs[i])
	Nsub = 1

	## Check Y12 spectra
	if sext.slit_width==3:
		## Ns
		Nsub = 2
	elif sext.slit_width==2:
		## Nh
		Nsub = 6
	
	## MC add pointing unc
	for j in range(Nmc+1):
		if j==0:
			sext.spec_build(fits_irc[i], Ny=24, Nsub=Nsub)
		else:
			sext.spec_build(fits_irc[i]+'_'+str(j), Ny=24, Nsub=Nsub)
'''
##---------------------------
##        Combine slits
##---------------------------
swp = iswarp(fits_irc, \
	# center='9:55:52,69:40:45', \
	pixscale=1.5, tmpdir=path_build, \
	verbose=verbose)
'''
## Reprendre MC adding spec unc
##------------------------------
for j in trange(Nmc+1, #leave=False, 
	desc='<iswarp> IRC Combining [MC]'):
	if j==0:
		comb = swp.combine(fits_irc, \
			'wgt_avg', keepedge=True, \
			tmpdir=path_build+'MC_no/', \
			filOUT=path_out+src+'_IRC_0')
	else:
		fits_irc_mc = []
		for f in fits_irc:
			fits_irc_mc.append(f+'_'+str(j))
		comb = swp.combine(fits_irc_mc, \
			keepedge=True, uncpdf='splitnorm', \
			tmpdir=path_build+'MC_'+str(j)+'/', \
			filOUT=path_out+src+'_IRC_'+str(j))

## Cal unc
##---------
mcimage = []
for j in trange(Nmc+1, #leave=False, 
	desc='IRC Reading [MC]'):
	if j==0:
		hd0 = read_fits(path_out+src+'_IRC_0')
		header = hd0.header
		wvl = hd0.wave
	else:
		hd = read_fits(path_out+src+'_IRC_'+str(j))
		mcimage.append(hd.data)
if Nmc>1:
	mcimage = np.array(mcimage)
	unc = np.nanstd(mcimage, axis=0)
	write_fits(path_out+src+'_IRC_unc', header, unc, wvl)
'''
##---------------------------
##     Inter-Calibration
##---------------------------

## Reproject phot
##----------------
image_phot = swp.combine([fits_phot], \
	keepedge=True, tmpdir=path_tmp, \
	filOUT=path_cal+src+'_'+phot).image[0]
image_phot0 = swp.combine([fits_phot0], \
	keepedge=True, tmpdir=path_tmp, \
	filOUT=path_cal+src+'_'+phot0).image[0]

## Convert unit
##--------------
cdelt = get_pc(header=read_fits(fits_phot).header).cdelt
unit_fac = gmean(1.e-6/pix2sr(1., cdelt)) # Jy/pix to MJy/sr
# print(cdelt, unit_fac)
image_phot *= unit_fac
hdr = swp.refheader
hdr['BUNIT'] = 'MJy/sr'
write_fits(path_cal+src+'_'+phot, swp.refheader, image_phot)

## Convolve phot
##---------------
if phot=='IRAC1':
	phot_ker = path_ker+'Kernel_HiRes_IRAC_3.6_to_Gauss_06.0'
elif phot=='IRAC2':
	phot_ker = path_ker+'Kernel_HiRes_IRAC_4.5_to_Gauss_06.0'
list_phot = [path_cal+src+'_'+phot, path_cal+src+'_'+phot0]

for p in list_phot:
	conv = iconvolve(p, [phot_ker]*2, csv_ker, filOUT=p+'_conv')
	conv.do_conv(ipath=path_idl)

## Synthetic photometry
##----------------------

# calib = intercalib(path_out+src+'_IRC_0')
# wcen, Fsyn, Fsig = calib.synthetic_photometry([phot])
## Alternative
hd = read_fits(path_out+src+'_IRC_0')
w_spec = hd.wave
w_spec[:2] = .1 # exclude 1st 2 wvl
Fnu_spec = np.copy(hd.data)
Fnu_spec[:2] = 0
wcen, Fsyn, Fsig = intercalib().synthetic_photometry([phot], 
	w_spec=w_spec, Fnu_spec=Fnu_spec)

ma_zero = np.ma.array(Fsyn, mask=(Fsyn==0)).mask
Fsyn[ma_zero] = np.nan
write_fits(path_cal+src+'_'+phot+'_sp', swp.refheader, Fsyn[0])

##---------------------------
##            Plots
##---------------------------

## Inter-Calibration
##-------------------
pts_synt = read_fits(path_cal+src+'_'+phot+'_sp').data.reshape((-1,))
pts_phot = read_fits(path_cal+src+'_'+phot+'_conv').data.reshape((-1,))
pts_phot0 = read_fits(path_cal+src+'_'+phot0+'_conv').data.reshape((-1,))

fig, ax = plt.subplots(1, 2, sharex=True, figsize=(11,5))
ax[0].errorbar(pts_phot0, pts_phot, 
	c='r', fmt='o', ms=2.)
ax[0].set_xlabel(phot+' (SINGS)')
ax[0].set_ylabel(phot+' (DP)')

mask = ~np.ma.array(pts_synt, 
	mask=np.logical_or(np.isnan(pts_synt), np.isnan(pts_phot0))).mask
ax[1].errorbar(pts_phot0[mask], pts_synt[mask], \
	fmt='o', ms=2., label='observations')
ax[1].set_xlabel(phot+' (SINGS)')
ax[1].set_ylabel('IRC (synt phot)')

## Linear fit
popt, pcov = curve_fit(f_lin, pts_phot0[mask], pts_synt[mask])
calib_fac = 1. / popt[0]
print('Inter-Calibration ('+phot+') factor (full image) = ', calib_fac)
ax[1].errorbar(pts_phot0[mask], f_lin(pts_phot0[mask], *popt), \
	c='m', ls='-', label='calib fac = {:.4}'.format(calib_fac))
ax[1].legend(loc='upper left')

plt.subplots_adjust(wspace=.3)

## Spectra
##---------
rx, ry = 86, 59 # pixel coord
xmin, xmax = 2, 6 # wavelength range
ymin, ymax = 2e2, 2e3 # surface brightness range

data_irc = read_fits(path_out+src+'_IRC_0', path_out+src+'_IRC_unc')
x = data_irc.wave
y = data_irc.data
dy = data_irc.unc

fig2, ax2 = plt.subplots(figsize=(8,5))
ax2.errorbar(x[2:], y[2:,ry,rx], \
	c='b', ls='dotted', label='no calib')
ax2.errorbar(x[2:], y[2:,ry,rx]*calib_fac, dy[2:,ry,rx], \
	c='k', ecolor='r', capsize=1, label='{} calib'.format(phot))
ax2.errorbar(wcen, Fsyn[0][ry,rx], yerr=Fsig[0][0], \
	fmt='^', c='g', capsize=1, label='IRC (synt phot)')
ax2.errorbar(wcen, read_fits(path_cal+src+'_'+phot0+'_conv').data[ry,rx], \
	fmt='o', c='olive', label=phot0+' (SINGS)')
ax2.legend(loc='upper left')
ax2.set_title(src+'_{}_{}'.format(rx, ry))

ax2.set_xscale('symlog')
ax2.set_yscale('symlog')
ax2.set_xlim((xmin, xmax))
ax2.set_ylim((ymin, ymax))
ax2.set_xticks(np.arange(xmin,xmax,1), minor=False)
ax2.set_yticks(np.arange(ymin,ymax,ymin), minor=False)
ax2.set_xlabel('Wavelength (micron)')
ax2.set_ylabel('Surface brightness (MJy/sr)')
ax2.xaxis.set_major_formatter(ScalarFormatter())
ax2.yaxis.set_major_formatter(ScalarFormatter())
# ax2.xaxis.set_minor_formatter(LogFormatter())
# ax2.yaxis.set_minor_formatter(NullFormatter())
ax2.vlines(3.3, 0, 1e4, linestyles='-.', colors='grey')
ax2.vlines(3.4, 0, 1e4, linestyles='-.', colors='grey')

plt.show()

fig.savefig(path_cal+src+'_'+phot)
fig2.savefig(path_cal+src+'_'+phot+'_{}_{}'.format(rx, ry))
