#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""

Build stitched Spitzer/IRS spectral cube

"""

import logging, sys
# logging.basicConfig(stream=sys.stderr, level=logging.ERROR)
# print(logging.getLogger())
logging.disable(sys.maxsize) # disable IDL print
import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning) 

from tqdm import tqdm, trange

import os
import numpy as np
## laputan
from laputan.inout import fclean, fitsext, read_fits, write_fits
from laputan.imaging import iconvolve, iswarp, concatenate
from laputan.astrom import fixwcs

##---------------------------
##      Initialisation
##---------------------------
## Local
from buildinfo import ( src, Nmc, verbose,
                        path_tmp, out_irc,
                        path_idl, path_ker, path_conv,
                        fits_ker, csv_ker,
                        chnl, fits_irs, out_irs
)
Nch = len(chnl)
# Nmc = 2

##---------------------------
##       Combine obs
##---------------------------
# swp = iswarp(sum(fits_irs, []), pixscale=6., tmpdir=path_tmp)
refheader = fixwcs(out_irc+'_0'+fitsext).header
swp = iswarp(refheader=refheader,
# swp = iswarp(sum(fits_irs, []), refheader=refheader,
             tmpdir=path_tmp, verbose=verbose)

## Add MC unc
##------------
for i in trange(Nch, #leave=False,
                desc='<iswarp> IRS Combining ({} chnl)'.format(Nch)):
    for j in trange(Nmc+1, leave=False,
                    desc='<iswarp> IRS Combining [MC]'):
        if j==0:
            comb = swp.combine(fits_irs[i], combtype='wgt_avg',
                               keepedge=True, cropedge=False,
                               tmpdir=path_tmp+'MC_no/', 
                               filOUT=path_tmp+src+'_'+chnl[i])
        else:
            comb = swp.combine(fits_irs[i], combtype='wgt_avg', 
                               keepedge=True, cropedge=False, uncpdf='norm',
                               tmpdir=path_tmp+'MC_'+str(j)+'/',
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
                             kfile=fits_ker, klist=csv_ker, convdir=path_conv,
                             filOUT=path_tmp+src+'_'+chnl[i]+'_'+str(j)+'_conv')
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
            hd = read_fits(path_tmp+src+'_'+chnl[i]+'_'+str(j)+'_conv')
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
        # filex = [out_irc+'_0']
        filex.extend(files)
        concatenate(filex, path_tmp+src+'_ma', wsort=False)
        data0 = read_fits(path_tmp+src+'_ma').data
        ma = np.ma.array(data0, mask=np.isnan(data0))
        mask_any = ma.mask.any(axis=0)
        mask_all = ma.mask.all(axis=0)

        concatenate(files, out_irs+'_0', wsort=False)
        hd = read_fits(out_irs+'_0')
        data = hd.data
        header = hd.header
        wvl = hd.wave
        for k in range(len(wvl)):
            data[k][mask_all] = np.nan
        write_fits(out_irs+'_0', header, data, wvl)
    else:
        for i in range(Nch):
            files.append(path_tmp+src+'_'+chnl[i]+'_'+str(j))

        concatenate(files, out_irs+'_'+str(j), wsort=False)
        data = read_fits(out_irs+'_'+str(j)).data
        for k in range(len(wvl)):
            data[k][mask_all] = np.nan
        write_fits(out_irs+'_'+str(j), header, data, wvl)

## Cal unc (all)
##---------------
mcimage = []
for j in trange(Nmc+1,
                desc='IRS Reading [MC]'):
    if j==0:
        hd0 = read_fits(out_irs+'_0')
        header = hd0.header
        wvl = hd0.wave
    else:
        hd = read_fits(out_irs+'_'+str(j))
        mcimage.append(hd.data)
if Nmc>1:
    mcimage = np.array(mcimage)
    unc = np.nanstd(mcimage, axis=0)
    write_fits(out_irs+'_unc', header, unc, wvl)
