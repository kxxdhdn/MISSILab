#!/usr/bin/env python
# -*- coding: utf-8 -*-

import logging, sys
# logging.basicConfig(stream=sys.stderr, level=logging.ERROR)
# print(logging.getLogger())
logging.disable(sys.maxsize)
import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning) 

from tqdm import tqdm, trange

import os
import numpy as np
from scipy.optimize import curve_fit
from scipy.stats import gmean
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter, LogFormatter, NullFormatter
## astylo
from astylo.bio import read_fits, write_fits
from astylo.proc import iconvolve, iswarp, concatenate
from astylo.calib import intercalib
from astylo.lib import fclean, f_lin, f_lin0
from astylo.astrolib import pix2sr, get_pc

## Local
from param import (
    src, Nmc, path_irs, path_ker, path_idl, 
    fits_irs, chnl, Nch, fits_ker, csv_ker, 
    phot, phot0, path_phot, path_cal, fits_phot, fits_phot0, 
    path_out, path_tmp, path_slices, path_tests, verbose, 
)

##---------------------------
##      Initialisation
##---------------------------
Nmc = 2

##---------------------------
##       Combine obs
##---------------------------
swp = iswarp(sum(fits_irs, []), \
    # center='9:55:52,69:40:45', \
    pixscale=1.5, tmpdir=path_tmp, \
    verbose=False)

## Reprendre MC adding spec unc
##------------------------------
for i in trange(Nch, #leave=False, 
    desc='<iswarp> IRS Combining ({} chnl)'.format(Nch)):
    for j in trange(Nmc+1, leave=False, 
        desc='<iswarp> IRS Combining [MC]'):
        if j==0:
            comb = swp.combine(fits_irs[i], \
                'wgt_avg', keepedge=True, \
                tmpdir=path_tmp+'MC_no/', \
                filOUT=path_tmp+src+'_'+chnl[i])
        else:
            comb = swp.combine(fits_irs[i], 'wgt_avg', \
                keepedge=True, uncpdf='norm', \
                tmpdir=path_tmp+'MC_'+str(j)+'/', \
                filOUT=path_tmp+src+'_'+chnl[i]+'_'+str(j))

## PSF Convolution
##-----------------
for i in trange(Nch, #leave=False, 
    desc='<iconvolve> IRS Smoothing ({} chnl)'.format(Nch)):
    for j in trange(Nmc+1, leave=False, 
        desc='<iconvolve> IRS Smoothing [MC]'):
        if j==0:
            conv = iconvolve(path_tmp+src+'_'+chnl[i], \
                fits_ker, csv_ker, \
                filTMP=path_slices+src+'_'+chnl[i], \
                filOUT=path_tmp+src+'_'+chnl[i]+'_conv')
        else:
            conv = iconvolve(path_tmp+src+'_'+chnl[i]+'_'+str(j), \
                fits_ker, csv_ker, \
                filTMP=path_slices+src+'_'+chnl[i]+'_'+str(j), \
                filOUT=path_tmp+src+'_'+chnl[i]+'_conv_'+str(j))
        conv.do_conv(ipath=path_idl)

## Cal unc (chnl)
##----------------
for i in trange(Nch, #leave=False, 
    desc='IRS Cal unc ({} chnl)'.format(Nch)):
    mcimage = []
    for j in trange(Nmc+1, leave=False, 
        desc='IRS Reading [MC]'):
        if j==0:
            hd0 = read_fits(path_tmp+src+'_'+chnl[i]+'_conv')
            header = hd0.header
            wvl = hd0.wave
        else:
            hd = read_fits(path_tmp+src+'_'+chnl[i]+'_conv_'+str(j))
            mcimage.append(hd.data)
    if Nmc>1:
        mcimage = np.array(mcimage)
        unc = np.nanstd(mcimage, axis=0)
        write_fits(path_tmp+src+'_'+chnl[i]+'_conv_unc', header, unc, wvl)

##---------------------------
##     Concatenate IRS
##---------------------------

for j in trange(Nmc+1, #leave=False, 
    desc='IRS concatenation [MC]'):
    files = []
    func = []
    if j==0:
        for i in range(Nch):
            files.append(path_tmp+src+'_'+chnl[i])
        for f in files:
            func.append(f+'_unc')

        ## Make NaN mask
        filex = []
        # filex = [path_out+src+'_IRC_0']
        filex.extend(files)
        concatenate(filex, path_tmp+src+'_ma', sort_wave=True)
        data0 = read_fits(path_tmp+src+'_ma').data
        ma = np.ma.array(data0, mask=np.isnan(data0))
        mask_any = ma.mask.any(axis=0)

        concatenate(files, path_out+src+'_IRS_0', sort_wave=True)
        hd = read_fits(path_out+src+'_IRS_0')
        data = hd.data
        header = hd.header
        wvl = hd.wave
        for k in range(len(wvl)):
            data[k][mask_any] = np.nan
        write_fits(path_out+src+'_IRS_0', header, data, wvl)
    else:
        for i in range(Nch):
            files.append(path_tmp+src+'_'+chnl[i]+'_'+str(j))

        concatenate(files, path_out+src+'_IRS_'+str(j), sort_wave=True)
        data = read_fits(path_out+src+'_IRS_'+str(j)).data
        for k in range(len(wvl)):
            data[k][mask_any] = np.nan
        write_fits(path_out+src+'_IRS_'+str(j), header, data, wvl)

## Cal unc (all)
##---------------
mcimage = []
for j in trange(Nmc+1, #leave=False, 
    desc='IRS Reading [MC]'):
    if j==0:
        hd0 = read_fits(path_out+src+'_IRS_0')
        header = hd0.header
        wvl = hd0.wave
    else:
        hd = read_fits(path_out+src+'_IRS_'+str(j))
        mcimage.append(hd.data)
if Nmc>1:
    mcimage = np.array(mcimage)
    unc = np.nanstd(mcimage, axis=0)
    write_fits(path_out+src+'_IRS_unc', header, unc, wvl)

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
if phot=='IRAC3':
    phot_ker = path_ker+'Kernel_HiRes_IRAC_5.8_to_Gauss_06.0'
elif phot=='IRAC4':
    phot_ker = path_ker+'Kernel_HiRes_IRAC_8.0_to_Gauss_06.0'
list_phot = [path_cal+src+'_'+phot, path_cal+src+'_'+phot0]

for p in list_phot:
    conv = iconvolve(p, [phot_ker]*2, csv_ker, filOUT=p+'_conv')
    conv.do_conv(ipath=path_idl)

## Synthetic photometry
##----------------------

# calib = intercalib(path_out+src+'_IRS_0')
# wcen, Fsyn, Fsig = calib.synthetic_photometry([phot])
## Alternative
hd = read_fits(path_out+src+'_IRS_0')
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
ax[1].set_ylabel('IRS (synt phot)')

## Linear fit
popt, pcov = curve_fit(f_lin0, pts_phot0[mask], pts_synt[mask])
calib_fac = 1. / popt[0]
print('Inter-Calibration ('+phot+') factor (full image) = ', calib_fac)
ax[1].errorbar(pts_phot0[mask], f_lin0(pts_phot0[mask], *popt), \
    c='m', ls='-', label='calib fac = {:.4}'.format(calib_fac))
ax[1].legend(loc='upper left')

plt.subplots_adjust(wspace=.3)

## Spectra
##---------
rx, ry = 254, 217 # pixel coord
xmin, xmax = 5, 22 # wavelength range
ymin, ymax = 5e2, 3e4 # surface brightness range

data_irs = read_fits(path_out+src+'_IRS_0', path_out+src+'_IRS_unc')
x = data_irs.wave
y = data_irs.data
dy = data_irs.unc

fig2, ax2 = plt.subplots(figsize=(8,5))
ax2.errorbar(x[:], y[:,ry,rx], \
    c='b', ls='dotted', label='no calib')
ax2.errorbar(x[:], y[:,ry,rx]*calib_fac, dy[:,ry,rx], \
    c='k', ecolor='r', capsize=1, label='{} calib'.format(phot))
ax2.errorbar(wcen, Fsyn[0][ry,rx], yerr=Fsig[0][0], \
    fmt='^', c='g', capsize=1, label='IRS (synt phot)')
ax2.errorbar(wcen, read_fits(path_cal+src+'_'+phot0+'_conv').data[ry,rx], \
    fmt='o', c='olive', label=phot+' (SINGS)')
ax2.legend(loc='upper left')
ax2.set_title(src+'_{}_{}'.format(rx, ry))

ax2.set_xscale('symlog')
ax2.set_yscale('symlog')
ax2.set_xlim((xmin, xmax))
ax2.set_ylim((ymin, ymax))
ax2.set_xticks(np.arange(xmin,xmax,1), minor=False)
ax2.set_yticks(np.arange(ymin,ymax,ymin*10), minor=False)
ax2.set_xlabel('Wavelength (micron)')
ax2.set_ylabel('Surface brightness (MJy/sr)')
ax2.xaxis.set_major_formatter(ScalarFormatter())
ax2.yaxis.set_major_formatter(ScalarFormatter())
# ax2.xaxis.set_minor_formatter(LogFormatter())
# ax2.yaxis.set_minor_formatter(NullFormatter())
ax2.vlines(6.2, 0, 5e4, linestyles='-.', colors='grey')
ax2.vlines(7.7, 0, 5e4, linestyles='-.', colors='grey')
ax2.vlines(8.6, 0, 5e4, linestyles='-.', colors='grey')
ax2.vlines(11.3, 0, 5e4, linestyles='-.', colors='grey')
ax2.vlines(17, 0, 5e4, linestyles='-.', colors='grey')

plt.show()

fig.savefig(path_cal+src+'_'+phot)
fig2.savefig(path_cal+src+'_'+phot+'_{}_{}'.format(rx, ry))
