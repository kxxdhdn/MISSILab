#!/usr/bin/env python
# -*- coding: utf-8 -*-

### Basic version

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
from astylo.iolib import fclean, read_fits, write_fits
from astylo.ipro import iconvolve, iswarp, concatenate
from astylo.calib import intercalib
from astylo.mlib import f_lin, f_lin0
from astylo.alib import pix2sr, get_pc

##---------------------------
##      Initialisation
##---------------------------
## Local
from param_irs_M83 import (
    verbose, src, Nmc, path_idl, path_irs, path_ker,
    chnl, Nch, fits_irs, fits_ker, csv_ker,
    path_tmp, path_conv, path_out,  path_tests, 
    phot, path_phot, path_cal, 
)

Nmc = 2

##---------------------------
##       Combine obs
##---------------------------
swp = iswarp(sum(fits_irs, []), pixscale=5., tmpdir=path_tmp)
'''
## Add MC unc
##------------
for i in trange(Nch, #leave=False,
                desc='<iswarp> IRS Combining ({} chnl)'.format(Nch)):
    for j in trange(Nmc+1, leave=False,
                    desc='<iswarp> IRS Combining [MC]'):
        if j==0:
            comb = swp.combine(fits_irs[i], \
                combtype='wgt_avg', keepedge=True, \
                tmpdir=path_tmp+'MC_no/', \
                filOUT=path_tmp+src+'_'+chnl[i])
        else:
            comb = swp.combine(fits_irs[i], 'wgt_avg', \
                keepedge=True, uncpdf='norm', \
                tmpdir=path_tmp+'MC_'+str(j)+'/', \
                filOUT=path_tmp+src+'_'+chnl[i]+'_'+str(j))

## PSF Convolution
##-----------------
for i in trange(Nch,
                desc='<iconvolve> IRS Smoothing ({} chnl)'.format(Nch)):
    for j in trange(Nmc+1, leave=False,
                    desc='<iconvolve> IRS Smoothing [MC]'):
        if j==0:
            conv = iconvolve(path_tmp+src+'_'+chnl[i],
                             kfile=fits_ker, klist=csv_ker, convdir=path_conv,
                             filOUT=path_tmp+src+'_'+chnl[i]+'_conv')
        else:
            conv = iconvolve(path_tmp+src+'_'+chnl[i]+'_'+str(j),
                             fits_ker, csv_ker, convdir=path_conv,
                             filOUT=path_tmp+src+'_'+chnl[i]+str(j)+'_conv')
        conv.do_conv(idldir=path_idl)

## Cal unc (chnl)
##----------------
for i in trange(Nch,
                desc='IRS Cal unc ({} chnl)'.format(Nch)):
    mcimage = []
    for j in trange(Nmc+1, leave=False,
                    desc='IRS Reading [MC]'):
        if j==0:
            hd0 = read_fits(path_tmp+src+'_'+chnl[i]+'_conv')
            header = hd0.header
            wvl = hd0.wave
        else:
            hd = read_fits(path_tmp+src+'_'+chnl[i]+str(j)+'_conv')
            mcimage.append(hd.data)
    if Nmc>1:
        mcimage = np.array(mcimage)
        unc = np.nanstd(mcimage, axis=0)
        write_fits(path_tmp+src+'_'+chnl[i]+'_conv_unc', header, unc, wvl)

##---------------------------
##     Concatenate IRS
##---------------------------
for j in trange(Nmc+1,
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
for j in trange(Nmc+1,
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
'''



##---------------------------
##     Inter-Calibration
##---------------------------
raw_phot1 = path_phot+src+'_'+phot+'_DP'
# raw_phot2 = path_phot+src+'_'+phot+'_SINGS'
fits_phot1 = path_cal+src+'_'+phot+'_DP'
# fits_phot2 = path_cal+src+'_'+phot+'_SINGS'
## Reproject phot
##----------------
image_phot1 = swp.combine(raw_phot1, tmpdir=path_tmp,
                          filOUT=fits_phot1).image
# image_phot2 = swp.combine(raw_phot2, tmpdir=path_tmp,
#                           filOUT=fits_phot2).image

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
if phot=='IRAC3':
    phot_ker = path_ker+'Kernel_HiRes_IRAC_5.8_to_Gauss_06.0'
elif phot=='IRAC4':
    phot_ker = path_ker+'Kernel_HiRes_IRAC_8.0_to_Gauss_06.0'
# list_phot = [fits_phot1, fits_phot2]
list_phot = [fits_phot1]

for p in list_phot:
    conv = iconvolve(p, [phot_ker]*2, csv_ker, filOUT=p+'_conv')
    conv.do_conv(path_idl)

## Synthetic photometry
##----------------------
calib = intercalib(path_out+src+'_IRS_0')
sp = calib.synthetic_photometry(phot)

write_fits(path_cal+src+'_'+phot+'_sp', swp.refheader, sp.Fnu_filt)

##---------------------------
##           Plots
##---------------------------

## Inter-Calibration
##-------------------
pts_synt = read_fits(path_cal+src+'_'+phot+'_sp').data.reshape((-1,))
pts_phot1 = read_fits(fits_phot1+'_conv').data.reshape((-1,))
# pts_phot2 = read_fits(fits_phot2+'_conv').data.reshape((-1,))

nrows, ncols = 1, 1
fig, axes = plt.subplots(nrows, ncols, sharex=False, figsize=(5,5))
if nrows==1 and ncols==1:
    ax = axes
else:
    ax = axes[0]
    
## src_phot left fig
mask0 = ~np.ma.array(pts_synt,
                     mask=np.logical_or(np.isnan(pts_synt),
                                        np.isnan(pts_phot1),)).mask
ax.errorbar(pts_phot1[mask0], pts_synt[mask0],
               fmt='o', ms=2., label='observations')
ax.set_xlabel(phot+' (DP)')
ax.set_ylabel('IRS (synt phot)')
ax.set_xlim(0,1.1*max(pts_phot1[mask0]))
ax.set_ylim(0,1.1*max(pts_synt[mask0]))
## Linear fit
popt, pcov = curve_fit(f_lin0, pts_phot1[mask0], pts_synt[mask0])
calib_fac = 1. / popt[0]
print('Inter-Calibration ('+phot+') factor (full image) = ', calib_fac)
ax.errorbar(pts_phot1[mask0], f_lin0(pts_phot1[mask0], *popt),
               c='m', ls='-', label='calib fac = {:.4}'.format(calib_fac))
ax.legend(loc='upper left')

## src_phot right fig
# Fnu_max = 600 # (MJy/sr), center region might have saturation
# mask1 = ~np.ma.array(pts_phot2,
#                      mask=pts_phot2>Fnu_max).mask
# axes[1].errorbar(pts_phot2[mask1], pts_phot1[mask1], c='r', fmt='o', ms=2.)
# axes[1].set_xlabel(phot+' (SINGS)')
# axes[1].set_ylabel(phot+' (DP)')
# axes[1].set_xlim((0,600))
# axes[1].set_ylim((0,600))

plt.subplots_adjust(wspace=.3)

## Spectra
##---------
rx, ry = 113, 85 # pixel coord
xmin, xmax = 5, 22 # wavelength range
ymin, ymax = 5e2, 3e4 # surface brightness range

data_irs = read_fits(path_out+src+'_IRS_0', path_out+src+'_IRS_unc')
x = data_irs.wave
y = data_irs.data
dy = data_irs.unc

fig2, ax2 = plt.subplots(figsize=(8,5))
ax2.errorbar(x[:], y[:,ry,rx],
             c='b', ls='dotted', label='no calib')
ax2.errorbar(x[:], y[:,ry,rx]*calib_fac, dy[:,ry,rx],
             c='k', ecolor='r', capsize=1, label='{} calib'.format(phot))
ax2.errorbar(sp.wcen, sp.Fnu_filt[ry,rx], yerr=sp.smat,
             fmt='^', c='g', capsize=1, label='IRS (synt phot)')
ax2.errorbar(sp.wcen, read_fits(fits_phot1+'_conv').data[ry,rx],
             fmt='o', c='olive', label=phot+' (DP)')
ax2.legend(loc='upper left')
ax2.set_title(src+'_{}_{}'.format(rx, ry))

ax2.set_xscale('symlog')
ax2.set_yscale('symlog')
ax2.set_xlim((xmin, xmax))
# ax2.set_ylim((ymin, ymax))
ax2.set_xticks(np.arange(xmin,xmax,1), minor=False)
ax2.set_yticks(np.arange(ymin,ymax,ymin*10), minor=False)
ax2.set_xlabel('Wavelength (micron)')
ax2.set_ylabel('Surface brightness (MJy/sr)')
ax2.xaxis.set_major_formatter(ScalarFormatter())
ax2.yaxis.set_major_formatter(ScalarFormatter())
# ax2.xaxis.set_minor_formatter(LogFormatter())
# ax2.yaxis.set_minor_formatter(NullFormatter())
ax2.vlines(6.2, .5*min(y[:,ry,rx]), 2*max(y[:,ry,rx]),
           linestyles='-.', colors='grey')
ax2.vlines(7.7, .5*min(y[:,ry,rx]), 2*max(y[:,ry,rx]),
           linestyles='-.', colors='grey')
ax2.vlines(8.6, .5*min(y[:,ry,rx]), 2*max(y[:,ry,rx]),
           linestyles='-.', colors='grey')
ax2.vlines(11.3, .5*min(y[:,ry,rx]), 2*max(y[:,ry,rx]),
           linestyles='-.', colors='grey')
ax2.vlines(17, .5*min(y[:,ry,rx]), 2*max(y[:,ry,rx]),
           linestyles='-.', colors='grey')

# plt.show()

fig.savefig(path_cal+src+'_'+phot)
fig2.savefig(path_cal+src+'_'+phot+'_{}_{}'.format(rx, ry))

##---------------------------
##           Clean
##---------------------------
# if input('Clean tmp (y/n): ')=='y':
#     fclean(path_tmp)
