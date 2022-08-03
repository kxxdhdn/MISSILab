#!/usr/bin/env python
# -*- coding: utf-8 -*-

# import logging, sys
# logging.disable(sys.maxsize)
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
from astylo.ipro import hextract, hswarp, iconvolve, iswarp, concatenate
from astylo.calib import intercalib
from astylo.mathlib import f_lin, f_lin0
from astylo.arrlib import closest
from astylo.astrolib import pix2sr, get_pc, fixwcs
from astylo.plotlib import plot2D, plot2D_m

## Local
from param import (
    src, Nmc, path_idl, path_ker, fits_ker, csv_ker, 
    phot, phot0, path_phot, path_cal, fits_phot, fits_phot0, 
    path_out, path_tmp, path_slices, path_tests, verbose, 
)

##---------------------------
##       Initialisation
##---------------------------
Nmc = 100
'''
##---------------------------
##     Rescale pixel size
##---------------------------
ds0 = read_fits(path_out+src+'_IRC_0')
oldimage = ds0.data
wvl = ds0.wave
oldheader = fixwcs(path_out+src+'_IRC_0').header
write_fits(path_tmp+src+'_IRC_fp', oldheader, oldimage[0])
'''
swp = iswarp([path_tmp+src+'_IRC_fp'], \
    pixscale=6., tmpdir=path_tmp, \
    verbose=False)
'''
swp.footprint(path_tmp+'footprint')

for j in trange(Nmc+1, #leave=False, 
    desc='Rescaling pixel size [MC]'):
    ## IRC
    swp.combine([path_out+src+'_IRC_'+str(j)], \
        keepedge=False, tmpdir=path_slices, \
        filOUT=path_tmp+src+'_IRC_'+str(j))
    ## IRS
    swp.combine([path_out+src+'_IRS_'+str(j)], \
        keepedge=False, tmpdir=path_slices, \
        filOUT=path_tmp+src+'_IRS_'+str(j))

## Cal unc (instr)
##-----------------
mcimage = []
for j in trange(Nmc+1, #leave=False, 
    desc='IRC Reading [MC]'):
    if j==0:
        hd0 = read_fits(path_tmp+src+'_IRC_0')
        header = hd0.header
        wvl = hd0.wave
    else:
        hd = read_fits(path_tmp+src+'_IRC_'+str(j))
        mcimage.append(hd.data)
if Nmc>1:
    mcimage = np.array(mcimage)
    unc = np.nanstd(mcimage, axis=0)
    write_fits(path_tmp+src+'_IRC_unc', header, unc, wvl)

mcimage = []
for j in trange(Nmc+1, #leave=False, 
    desc='IRS Reading [MC]'):
    if j==0:
        hd0 = read_fits(path_tmp+src+'_IRS_0')
        header = hd0.header
        wvl = hd0.wave
    else:
        hd = read_fits(path_tmp+src+'_IRS_'+str(j))
        mcimage.append(hd.data)
if Nmc>1:
    mcimage = np.array(mcimage)
    unc = np.nanstd(mcimage, axis=0)
    write_fits(path_tmp+src+'_IRS_unc', header, unc, wvl)

##---------------------------
##       Concatenate ALL
##---------------------------
for j in trange(Nmc+1, #leave=False, 
    desc='ALL concatenation [MC]'):
    files = [path_tmp+src+'_IRC_'+str(j), path_tmp+src+'_IRS_'+str(j)]
    if j==0:
        func = [path_tmp+src+'_IRC_unc', path_tmp+src+'_IRS_unc']

        concatenate(files, path_tmp+src+'_ma', sort_wave=True)
        data0 = read_fits(path_tmp+src+'_ma').data
        ma = np.ma.array(data0, mask=np.isnan(data0))
        mask_any = ma.mask.any(axis=0)

        concatenate(files, path_out+src+'_0', sort_wave=True)
        ds = read_fits(path_out+src+'_0')
        data = ds.data
        header = ds.header
        wvl = ds.wave
        for k in range(len(wvl)):
            data[k][mask_any] = np.nan
        write_fits(path_out+src+'_0', header, data, wvl)
    else:
        concatenate(files, path_out+src+'_'+str(j), sort_wave=True)
        data = read_fits(path_out+src+'_'+str(j)).data
        for k in range(len(wvl)):
            data[k][mask_any] = np.nan
        write_fits(path_out+src+'_'+str(j), header, data, wvl)

## Cal unc (ALL)
##---------------
mcimage = []
for j in trange(Nmc+1, #leave=False, 
    desc='ALL Reading [MC]'):
    if j==0:
        hd0 = read_fits(path_out+src+'_0')
        header = hd0.header
        wvl = hd0.wave
    else:
        hd = read_fits(path_out+src+'_'+str(j))
        mcimage.append(hd.data)
if Nmc>1:
    mcimage = np.array(mcimage)
    unc = np.nanstd(mcimage, axis=0)
    write_fits(path_out+src+'_unc', header, unc, wvl)
'''

##---------------------------
##       Plot spectra
##---------------------------
phot = ['IRAC1','IRAC2','IRAC4','MIPS1']

## Synthetic photometry
##----------------------
calib = intercalib(path_out+src+'_0')
sp = []
for i in range(4):
    raw_phot2 = path_phot+src+'_'+phot[i]+'_SINGS'
    fits_phot2 = path_cal+src+'_'+phot[i]+'_SINGS'
    ## Reproject phot
    ##----------------
    image_phot2 = swp.combine(raw_phot2, tmpdir=path_tmp,
                              filOUT=fits_phot2).image
    
    ## Convolve phot
    ##---------------
    if phot[i]=='IRAC1':
        phot_ker = path_ker+'Kernel_HiRes_IRAC_3.6_to_Gauss_06.0'
    elif phot[i]=='IRAC2':
        phot_ker = path_ker+'Kernel_HiRes_IRAC_4.5_to_Gauss_06.0'
    elif phot[i]=='IRAC3':
        phot_ker = path_ker+'Kernel_HiRes_IRAC_5.8_to_Gauss_06.0'
    elif phot[i]=='IRAC4':
        phot_ker = path_ker+'Kernel_HiRes_IRAC_8.0_to_Gauss_06.0'
    elif phot[i]=='MIPS1':
        phot_ker = path_ker+'Kernel_HiRes_MIPS_24_to_Gauss_06.0'
    list_phot = [fits_phot2]
    
    for p in list_phot:
        conv = iconvolve(p, [phot_ker]*2, csv_ker, filOUT=p+'_conv')
        conv.do_conv(idldir=path_idl)
    
    sp.append(calib.synthetic_photometry(phot[i]))

    write_fits(path_cal+src+'_'+phot[i]+'_sp', swp.refheader, sp[i].Fnu_filt)

## Read data
syn_map = []
pho_map = []
syn = []
pho = []
calib_fac = []
calib_off = []
for i in range(4):
    syn_map.append(read_fits(path_cal+src+'_'+phot[i]+'_sp').data)
    pho_map.append(read_fits(path_cal+src+'_'+phot[i]+'_SINGS_conv').data)
    syn.append(syn_map[i].reshape((-1,)))
    pho.append(pho_map[i].reshape((-1,)))

nrows, ncols = 2, 2
fig, axes = plt.subplots(nrows, ncols, sharex=False, figsize=(8,8))
for i in range(4):
    if i<2:
        mask0 = ~np.ma.array(syn[i],
                         mask=np.logical_or(
                           np.logical_or(
                             np.logical_or(
                               np.isnan(syn[i]),np.isnan(pho[i])),
                                   syn[i]>6.e1),
                                 pho[i]>6.e1)).mask
    else:
        mask0 = ~np.ma.array(syn[i],
                         mask=np.logical_or(
                           np.logical_or(
                             np.logical_or(
                               np.isnan(syn[i]),np.isnan(pho[i])),
                                   syn[i]>1.e3),
                                 pho[i]>1.e3)).mask
    px,py = int(i/2), i%2
    axes[px,py].errorbar(pho[i][mask0], syn[i][mask0],
                     fmt='o', ms=2., label=src+' observations')
    axes[px,py].set_xlabel(phot[i]+' (SINGS)')
    axes[px,py].set_ylabel('Synthetic photometry')
    axes[px,py].set_xlim(0,1.1*max(syn[i][mask0]))
    axes[px,py].set_ylim(0,1.1*max(syn[i][mask0]))

    ## Linear fit
    popt, pcov = curve_fit(f_lin0, pho[i][mask0], syn[i][mask0])
    calib_fac.append(1. / popt[0])
    # calib_off.append(-popt[1])
    print('Inter-Calibration ('+phot[i]+') factor (full image) = ', calib_fac[i])
    # print('Inter-Calibration ('+phot[i]+') offset (full image) = ', calib_off[i])
    axes[px,py].errorbar(pho[i][mask0], f_lin0(pho[i][mask0], *popt),
                   c='m', ls='-', label='calib fac = {:.4}'.format(calib_fac[i]))
    axes[px,py].legend(loc='upper left')

fig.savefig(path_cal+src+'_intercalib')


## Spectra
##---------
# x = [6,12,14,18,18,18,19,19,20,20,22,49,49,50,51,52,68,]
# y = [15,20,21,28,22,23,22,23,20,21,21,43,47,46,48,45,118,]

rx, ry = 22,21
ymin, ymax = 5, 1e4 # surface brightness range
# rx, ry = 52,45
# ymin, ymax = -1, 1e1 # surface brightness range
xmin, xmax = 2., 40. # wavelength range

## Read data
ds = read_fits(path_out+src+'_0', path_out+src+'_unc')
data = ds.data
unc = ds.unc
# data[:259,:,:] = data[:259,:,:]*data[259,:,:]/data[258,:,:]
wvl = ds.wave
header = ds.header

iwdel = []
iwdel.append(closest(wvl,5.))
iwdel.append(closest(wvl,5.))
iwdel.append(closest(wvl,5.))
wvl = np.delete(wvl,iwdel)
data = np.delete(data,iwdel,axis=0)
unc = np.delete(unc,iwdel,axis=0)
# print(wvl[257:260])

wrange = [ (2.50, 5.00), # irc - IRAC1 (3.08-4.01)
           (5.21, 14.28), # sl2+sl1 - IRAC4 (6.15-10.50)
           (14.29, 38.00), ] # ll2+ll1 - MIPS1 (18.01-32.21)
iwran = []
for i in range(3):
    iwmin = closest(wvl,wrange[i][0])
    iwmax = closest(wvl,wrange[i][1])+1
    # print(iwmin,iwmax)
    iwran.append((iwmin,iwmax))

data_new = data.copy()
for i in range(3):
    if i==0:
        data_new[iwran[i][0]:iwran[i][1],:,:] = data[iwran[i][0]:iwran[i][1],:,:]*calib_fac[i]
    else:
        data_new[iwran[i][0]:iwran[i][1],:,:] = data[iwran[i][0]:iwran[i][1],:,:]*calib_fac[i+1]


write_fits(path_out+src+'_calib', header, data_new[:,21:23,22:26], wvl)
write_fits(path_out+src+'_calib_unc', header, unc[:,21:23,22:26], wvl)

fig2, ax2 = plt.subplots(figsize=(8,5))
ax2.errorbar(wvl, data[:,ry,rx], \
    c='b', ls='dotted', label='no calib')
ax2.errorbar(wvl, data_new[:,ry,rx], unc[:,ry,rx], \
    c='k', ecolor='r', capsize=1, label='calib')
ax2.errorbar(sp[0].wcen, syn_map[0][ry,rx], yerr=sp[0].smat,
             fmt='^', c='g', capsize=1, label='IRC (synt phot)')
ax2.errorbar(sp[0].wcen, pho_map[0][ry,rx],
             fmt='o', c='olive', label=phot[0]+' (SINGS)')
ax2.errorbar(sp[2].wcen, syn_map[2][ry,rx], yerr=sp[2].smat,
             fmt='^', c='g', capsize=1, label='IRS (synt phot)')
ax2.errorbar(sp[2].wcen, pho_map[2][ry,rx],
             fmt='o', c='olive', label=phot[2]+' (SINGS)')
ax2.errorbar(sp[3].wcen, syn_map[3][ry,rx], yerr=sp[3].smat,
             fmt='^', c='g', capsize=1, label='IRS (synt phot)')
ax2.errorbar(sp[3].wcen, pho_map[3][ry,rx],
             fmt='o', c='olive', label=phot[3]+' (SINGS)')
ax2.legend(loc='upper left')
ax2.set_title(src+'_{}_{}'.format(rx+1, ry+1))

ax2.set_xscale('symlog')
ax2.set_yscale('symlog')
ax2.set_xlim((xmin, xmax))
ax2.set_ylim((ymin, ymax))

# tics = [2,3,4,5,6,7,8,9,10,11,13,17,20,25,30,35,40]
tics = [2,3,4,5,6,7,8,9,10,11,12,13,14,15,17,20]
ax2.set_xticks(tics, minor=False)
# ax2.set_xticks(np.arange(xmin,xmax,1), minor=False)
# ax2.set_yticks(np.arange(ymin,ymax,ymin), minor=False)
ax2.set_xlabel('Wavelength (micron)')
ax2.set_ylabel('Surface brightness (MJy/sr)')
ax2.xaxis.set_major_formatter(ScalarFormatter())
ax2.yaxis.set_major_formatter(ScalarFormatter())
# ax2.xaxis.set_minor_formatter(LogFormatter())
# ax2.yaxis.set_minor_formatter(NullFormatter())
ax2.vlines(5., 0, 1e4, linestyles='-.', colors='grey')
ax2.vlines(5.21, 0, 1e4, linestyles='-.', colors='grey')
ax2.vlines(7.57, 0, 1e4, linestyles='-.', colors='grey')
ax2.vlines(14.29, 0, 1e4, linestyles='-.', colors='grey')
ax2.vlines(20.67, 0, 1e4, linestyles='-.', colors='grey')

fig2.savefig(path_cal+src+'_{}_{}'.format(rx+1, ry+1))



# plt.show()
