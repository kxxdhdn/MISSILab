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
from astylo.iolib import fclean, read_fits, write_fits
from astylo.ipro import iconvolve, iswarp, sextract
from astylo.calib import intercalib
from astylo.mlib import f_lin, f_lin0
from astylo.alib import pix2sr, get_pc

##---------------------------
##       Initialisation
##---------------------------
## Local
from param_irc import (
    verbose, src, Nmc, path_idl, path_irc, path_ker,
    parobs, fits_irc, csv_ker, 
    path_tmp, path_out, path_build, path_tests,
    phot, path_phot, path_cal, 
)
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
    pixscale=6, tmpdir=path_build, \
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
raw_phot1 = path_phot+src+'_'+phot+'_DP'
raw_phot2 = path_phot+src+'_'+phot+'_SINGS'
fits_phot1 = path_cal+src+'_'+phot+'_DP'
fits_phot2 = path_cal+src+'_'+phot+'_SINGS'
## Reproject phot
##----------------
image_phot1 = swp.combine(raw_phot1, tmpdir=path_tmp,
                          filOUT=fits_phot1).image
image_phot2 = swp.combine(raw_phot2, tmpdir=path_tmp,
                          filOUT=fits_phot2).image

## Convert phot (DP) unit
##------------------------
cdelt = get_pc(header=read_fits(raw_phot1).header).cdelt
unit_fac = gmean(1.e-6/pix2sr(1., cdelt)) # Jy/pix to MJy/sr
# print(cdelt, unit_fac)
image_phot1 *= unit_fac
hdr = swp.refheader
hdr['BUNIT'] = 'MJy/sr'
write_fits(fits_phot1, swp.refheader, image_phot1)

## Convolve phot
##---------------
if phot=='IRAC1':
    phot_ker = path_ker+'Kernel_HiRes_IRAC_3.6_to_Gauss_06.0'
elif phot=='IRAC2':
    phot_ker = path_ker+'Kernel_HiRes_IRAC_4.5_to_Gauss_06.0'
list_phot = [fits_phot1, fits_phot2]

for p in list_phot:
    conv = iconvolve(p, [phot_ker]*2, csv_ker, filOUT=p+'_conv')
    conv.do_conv(idldir=path_idl)

## Synthetic photometry
##----------------------
calib = intercalib(path_out+src+'_IRC_0')
sp = calib.synthetic_photometry(phot)

write_fits(path_cal+src+'_'+phot+'_sp', swp.refheader, sp.Fnu_filt)

##---------------------------
##            Plots
##---------------------------

## Inter-Calibration
##-------------------
pts_synt = read_fits(path_cal+src+'_'+phot+'_sp').data.reshape((-1,))
pts_phot1 = read_fits(fits_phot1+'_conv').data.reshape((-1,))
pts_phot2 = read_fits(fits_phot2+'_conv').data.reshape((-1,))

nrows, ncols = 1, 2
fig, axes = plt.subplots(nrows, ncols, sharex=False, figsize=(11,5))
if nrows==1 and ncols==1:
    ax = axes
else:
    ax = axes[0]
    
## src_phot left fig
mask0 = ~np.ma.array(pts_synt,
                     mask=np.logical_or(np.isnan(pts_synt),
                                        np.isnan(pts_phot2),)).mask
ax.errorbar(pts_phot2[mask0], pts_synt[mask0],
               fmt='o', ms=2., label='observations')
ax.set_xlabel(phot+' (SINGS)')
ax.set_ylabel('IRS (synt phot)')
ax.set_xlim(0,1.1*max(pts_phot2[mask0]))
ax.set_ylim(0,1.1*max(pts_synt[mask0]))
## Linear fit
popt, pcov = curve_fit(f_lin0, pts_phot2[mask0], pts_synt[mask0])
calib_fac = 1. / popt[0]
print('Inter-Calibration ('+phot+') factor (full image) = ', calib_fac)
ax.errorbar(pts_phot2[mask0], f_lin0(pts_phot2[mask0], *popt),
               c='m', ls='-', label='calib fac = {:.4}'.format(calib_fac))
ax.legend(loc='upper left')

## src_phot right fig
Fnu_max = 600 # (MJy/sr), center region might have saturation
mask1 = ~np.ma.array(pts_phot2,
                     mask=pts_phot2>Fnu_max).mask
axes[1].errorbar(pts_phot2[mask1], pts_phot1[mask1], c='r', fmt='o', ms=2.)
axes[1].set_xlabel(phot+' (SINGS)')
axes[1].set_ylabel(phot+' (DP)')
axes[1].set_xlim((0,600))
axes[1].set_ylim((0,600))

plt.subplots_adjust(wspace=.3)

## Spectra
##---------
rx, ry = 23, 16
ymin, ymax = 1e1, 2e2
# rx, ry = 86, 59 # pixel coord
xmin, xmax = 2, 6 # wavelength range
# ymin, ymax = 2e2, 2e3 # surface brightness range

data_irc = read_fits(path_out+src+'_IRC_0', path_out+src+'_IRC_unc')
x = data_irc.wave
y = data_irc.data
dy = data_irc.unc

fig2, ax2 = plt.subplots(figsize=(8,5))
ax2.errorbar(x[2:], y[2:,ry,rx], \
    c='b', ls='dotted', label='no calib')
ax2.errorbar(x[2:], y[2:,ry,rx]*calib_fac, dy[2:,ry,rx], \
    c='k', ecolor='r', capsize=1, label='{} calib'.format(phot))
ax2.errorbar(sp.wcen, sp.Fnu_filt[ry,rx], yerr=sp.smat,
             fmt='^', c='g', capsize=1, label='IRS (synt phot)')
ax2.errorbar(sp.wcen, read_fits(fits_phot2+'_conv').data[ry,rx],
             fmt='o', c='olive', label=phot+' (SINGS)')
ax2.legend(loc='upper left')
ax2.set_title(src+'_{}_{}'.format(rx+1, ry+1))

# ax2.set_xscale('symlog')
# ax2.set_yscale('symlog')
ax2.set_xlim((xmin, xmax))
ax2.set_ylim((ymin, ymax))
# ax2.set_xticks(np.arange(xmin,xmax,1), minor=False)
# ax2.set_yticks(np.arange(ymin,ymax,ymin), minor=False)
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
fig2.savefig(path_cal+src+'_'+phot+'_{}_{}'.format(rx+1, ry+1))
