#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""

Build AKARI/IRC spectral cube

"""

# import logging, sys
# logging.disable(sys.maxsize)
import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning) 
warnings.filterwarnings("ignore", message="Skipping SYSTEM_VARIABLE record")

from tqdm import tqdm, trange

import os
import numpy as np
## laputan
from laputan.inout import fclean, fitsext, read_fits, write_fits
from laputan.imaging import iconvolve, iswarp, sextract, icrop, irebin
from laputan.astrom import fixwcs

##---------------------------
##       Initialisation
##---------------------------
## Local
from buildinfo import ( src, Nmc, verbose,
                        path_irc, path_build,
                        parobs, fits_irc, out_irc
)
Nobs = len(parobs)
# Nmc = 2

##---------------------------
##        Build slits
##---------------------------

for i in trange(Nobs, #leave=False,
                desc='<sextract> IRC slit building'):
    sext = sextract(path_irc, parobs[i])
    Nsub = 6

    ## Reproduce Y12 spectra
    # if fits_irc[i][-1]=='s':
    #     Nsub = 2 # Ns
    # if fits_irc[i][-1]=='h':
    #     Nsub = 6 # Nh
    
    ## MC add pointing unc
    for j in range(Nmc+1):
        if j==0:
            if not os.path.isfile(fits_irc[i]+fitsext):
                sext.spec_build(fits_irc[i],
                                Nx=4, Ny=24, Nsub=Nsub)
                irebin(fits_irc[i], filOUT=fits_irc[i], pixscale=5)
                irebin(fits_irc[i]+'_unc', filOUT=fits_irc[i]+'_unc', pixscale=5)
                irebin(fits_irc[i]+'_unc_P', filOUT=fits_irc[i]+'_unc_P', pixscale=5)
                irebin(fits_irc[i]+'_unc_N', filOUT=fits_irc[i]+'_unc_N', pixscale=5)
        else:
            if not os.path.isfile(fits_irc[i]+'_'+str(j)+fitsext):
                sext.spec_build(fits_irc[i]+'_'+str(j),
                                Nx=4, Ny=24, Nsub=Nsub)
                irebin(fits_irc[i]+'_'+str(j), filOUT=fits_irc[i]+'_'+str(j),pixscale=5)
                irebin(fits_irc[i]+'_'+str(j)+'_unc', filOUT=fits_irc[i]+'_'+str(j)+'_unc', pixscale=5)
                irebin(fits_irc[i]+'_'+str(j)+'_unc_P', filOUT=fits_irc[i]+'_'+str(j)+'_unc_P', pixscale=5)
                irebin(fits_irc[i]+'_'+str(j)+'_unc_N', filOUT=fits_irc[i]+'_'+str(j)+'_unc_N', pixscale=5)

##---------------------------
##        Combine slits
##---------------------------
refheader = fixwcs(fits_irc[0]+fitsext).header
# swp = iswarp(refheader=refheader, \
swp = iswarp(fits_irc, refheader=refheader, 
             # center='9:55:52,69:40:45', pixscale=6,
             tmpdir=path_build, verbose=verbose)

## Reprendre MC adding spec unc
##------------------------------
for j in trange(Nmc+1, #leave=False, 
                desc='<iswarp> IRC Combining [MC]'):
    if j==0:
        comb = swp.combine(fits_irc, combtype='wgt_avg',
                           keepedge=True, cropedge=True,
                           tmpdir=path_build+'MC_no/',
                           filOUT=out_irc+'_0')
    else:
        fits_irc_mc = []
        for f in fits_irc:
            fits_irc_mc.append(f+'_'+str(j))
        comb = swp.combine(fits_irc_mc, combtype='wgt_avg',
                           keepedge=True, cropedge=True, uncpdf='splitnorm',
                           tmpdir=path_build+'MC_'+str(j)+'/',
                           filOUT=out_irc+'_'+str(j))

## Cal unc
##---------
mcimage = []
for j in trange(Nmc+1, #leave=False, 
                desc='IRC Reading [MC]'):
    if j==0:
        hd0 = read_fits(out_irc+'_0')
        header = hd0.header
        wvl = hd0.wave
    else:
        hd = read_fits(out_irc+'_'+str(j))
        mcimage.append(hd.data)
if Nmc>1:
    print('Calculating uncertainty cube...')
    mcimage = np.array(mcimage)
    unc = np.nanstd(mcimage, axis=0)
    write_fits(out_irc+'_unc', header, unc, wvl)
