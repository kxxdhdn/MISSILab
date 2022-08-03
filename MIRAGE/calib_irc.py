#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""

Inter-calibrate AKARI/IRC spectra with Spitzer/IRAC1

"""

# import logging, sys
# logging.disable(sys.maxsize)
import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning) 
warnings.filterwarnings("ignore", message="Skipping SYSTEM_VARIABLE record")

from tqdm import tqdm, trange

import os
import numpy as np
from scipy.optimize import curve_fit
## laputan
from laputan.inout import fclean, fitsext, read_fits, write_fits
from laputan.imaging import ( iconvolve, iswarp, icrop, iuncert,
                              respect, Jy_per_pix_to_MJy_per_sr )
from laputan.astrom import fixwcs, pix2sr
from laputan.calib import intercalib
from laputan.maths import f_lin, f_lin0, f_lin_p
from laputan.plots import pplot

##---------------------------
##       Preliminaries
##---------------------------
## Local
from buildinfo import ( src, Nmc, verbose,
                        path_phot, path_ker, path_cal,
                        path_idl, csv_ker, path_tmp,
                        out_irc, path_fig
)
# Nmc = 2

## Reprojection footprint
refheader = fixwcs(out_irc+'_0'+fitsext).header
swp = iswarp(refheader=refheader,
             tmpdir=path_tmp, verbose=verbose)

Nx = refheader['NAXIS1']
Ny = refheader['NAXIS2']

## Prepare photometry
##--------------------
phot = 'IRAC1'
phot_ker = path_ker+'Kernel_HiRes_IRAC_3.6_to_Gauss_06.0'

raw_phot1 = path_phot+src+'_'+phot+'_DP'
fits_phot1 = path_cal+src+'_'+phot+'_DP'

raw_phot2 = path_phot+src+'_'+phot+'_SINGS'
wgt_phot2 = path_phot+src+'_'+phot+'_SINGS_wgt'
fits_phot2 = path_cal+src+'_'+phot+'_SINGS'

fits_spec = path_cal+src+'_'+phot+'_IRC'

precalib = input("Run inter-calibration prepipeline (y/n)? ")
if precalib:
    ## DustPedia (phot1)
    ##===================
    
    ## Convert phot unit (Jy/pix -> MJy/sr)
    ## This step should be before iswarp reprojection (pixscale changed)
    Jy_per_pix_to_MJy_per_sr(raw_phot1, filOUT=fits_phot1)
    
    ## Reproject phot
    image_phot1 = swp.combine(fits_phot1,# combtype='wgt_avg',
                              cropedge=False, tmpdir=path_tmp,
                              filOUT=fits_phot1).image
    
    ## Convolve phot
    conv = iconvolve(fits_phot1, phot_ker, csv_ker, filOUT=fits_phot1)
    conv.do_conv(idldir=path_idl)
    
    ## DustPedia IRAC1 of M82 has no uncertainty map
    
    ## SINGS (phot2)
    ##===============
    
    ## Create uncertainty via weight map (suppose uniform contribution)
    ## The weight maps contain the information on the number of frames 
    ## that were used to create the science mosaics at each pixel (value= # frames x10); 
    ## the pixel size of the weight maps is the same as the science mosaics.
    ## -- SINGS v5 release doc
    bg = read_fits(raw_phot2).data[694:734,56:96]
    bg_wgt = read_fits(wgt_phot2).data[694:734,56:96] / 10.
    iuncert(raw_phot2, filOUT=raw_phot2+'_unc',
            filWGT=wgt_phot2, wfac=.1,
            BG_image=bg, BG_weight=bg_wgt)
    
    for j in trange(Nmc+1, leave=False,
                    desc='Generating SINGS uncertainty map [MC]'):
        if j==0:
            ## Reproject phot
            comb = swp.combine(raw_phot2, combtype='wgt_avg',
                               cropedge=False, tmpdir=path_tmp,
                               filOUT=fits_phot2).image
        
            ## Convolve phot
            conv = iconvolve(fits_phot2, phot_ker, csv_ker, filOUT=fits_phot2)
        else:
            comb = swp.combine(raw_phot2, combtype='wgt_avg',
                               cropedge=False, tmpdir=path_tmp, uncpdf='norm',
                               filOUT=path_tmp+src+'_'+phot+'_SINGS_'+str(j)).image
        
            conv = iconvolve(path_tmp+src+'_'+phot+'_SINGS_'+str(j),
                             kfile=phot_ker, klist=csv_ker,
                             filOUT=path_tmp+src+'_'+phot+'_SINGS_'+str(j))
        conv.do_conv(idldir=path_idl)
    
    mcimage = []
    for j in range(Nmc+1):
        if j==0:
            hdr = read_fits(fits_phot2).header
        else:
            mcimage.append(read_fits(path_tmp+src+'_'+phot+'_SINGS_'+str(j)).data)
    if Nmc>1:
        mcimage = np.array(mcimage)
        unc = np.nanstd(mcimage, axis=0)
        write_fits(fits_phot2+'_unc', hdr, unc)
    
    
    ## Synthetic photometry (spec)
    ##=============================
    mcimage = []
    for j in trange(Nmc+1,
                    desc='IRC sythetic photometry [MC]'):
        ic = intercalib(out_irc+'_'+str(j))
        sp = ic.synthetic_photometry(phot)
        if j==0:
            write_fits(fits_spec, ic.hdr, sp.Fnu_filt)
        else:
            ## Cal unc
            mcimage.append(sp.Fnu_filt)
    if Nmc>1:
        print('Calculating uncertainty cube...')
        mcimage = np.array(mcimage)
        unc = np.nanstd(mcimage, axis=0)
        write_fits(fits_spec+'_unc', ic.hdr, unc)

##---------------------------
##     Inter-Calibration
##---------------------------
data_phot1 = read_fits(fits_phot1).data
data_phot2 = read_fits(fits_phot2).data
unc_phot2 = read_fits(fits_phot2+'_unc').data
data_spec = read_fits(fits_spec).data
unc_spec = read_fits(fits_spec+'_unc').data

# clist = ['k', 'orange', 'g', 'c', 'b', 'm', 'pink']

Nz = 22
for z in range(Nz):
    if z==0:
        ## Zone 0 (all map)
        xmin, xmax = 0, Nx
        ymin, ymax = 0, Ny
        zname = 'all map'
    ## wind S
    elif z<3:
        ## Zone 1-2 (wind S 1-2)
        xmin, xmax = 53+z, 54+z
        ymin, ymax = 0, 12
        zname = 'wind S '+str(z)
    elif z<5:
        ## Zone 3-4 (wind S 3-4)
        xmin, xmax = 54+z, 55+z
        ymin, ymax = 0, 12
        zname = 'wind S '+str(z)
    elif z<9:
        ## Zone 5-8 (wind S 5-8)
        xmin, xmax = 55+z, 56+z
        ymin, ymax = 0, 12
        zname = 'wind S '+str(z)
    elif z<11:
        ## Zone 9-10 (center 1-2)
        xmin, xmax = 42+z, 43+z
        ymin, ymax = 12, 30
        zname = 'center '+str(z-8)
    elif z<13:
        ## Zone 11-12 (center 3-4)
        xmin, xmax = 44+z, 45+z
        ymin, ymax = 12, 30
        zname = 'center '+str(z-8)
    elif z==13:
        ## Zone 13 (center 5)
        xmin, xmax = 44+z, 48+z
        ymin, ymax = 12, 30
        zname = 'center '+str(z-8)
    elif z<18:
        ## Zone 14-17 (center 6-9)
        xmin, xmax = 46+z, 47+z
        ymin, ymax = 12, 30
        zname = 'center '+str(z-8)
    elif z==18:
        ## Zone 18 (inner wind N)
        xmin, xmax = 0, 53
        ymin, ymax = 30, 55
        zname = 'inner wind N'
    elif z==19:
        ## Zone 19 (wind N)
        xmin, xmax = 53, Nx
        ymin, ymax = 30, 55
        zname = 'wind N'
    elif z==20:
        ## Zone 20 (outer wind N)
        xmin, xmax = 0, Nx
        ymin, ymax = 55, 70
        zname = 'outer wind N'
    elif z==21:
        ## Zone 21 (cap N)
        xmin, xmax = 0, Nx
        ymin, ymax = 70, Ny
        zname = 'cap N'

        
    ## Plot
    ##------
    pix_phot1 = data_phot1[ymin:ymax,xmin:xmax].reshape((-1,))
    pix_phot2 = data_phot2[ymin:ymax,xmin:xmax].reshape((-1,))
    pix_spec = data_spec[ymin:ymax,xmin:xmax].reshape((-1,))
    pix_phot2_unc = unc_phot2[ymin:ymax,xmin:xmax].reshape((-1,))
    pix_spec_unc = unc_spec[ymin:ymax,xmin:xmax].reshape((-1,))
    
    # xgrid = np.arange(0,1e3,1)
    xgrid = np.logspace(-3,3,10000)

    ## NaNs mask
    mask = ~np.ma.array(pix_spec,
                        mask=np.logical_or(
                            np.isnan(pix_spec),np.isnan(pix_phot2))).mask

    ## S/N ratio
    # if z==0:
    #     print('S/N (IRC) = \n', pix_spec[mask]/pix_spec_unc[mask])
    #     print('S/N (SINGS) = \n', pix_phot2[mask]/pix_phot2_unc[mask])
    #     exit()
        
    if z<Nz-1:
        if z==0:
            ## SINGS - DP
            p0 = pplot(pix_phot1,  pix_phot2, yerr=pix_phot2_unc,
                       fmt='.', c='y', ec='r',elw=10,
                       xlog=1, ylog=1,
                       xlim=(0,1e3), ylim=(0,1e3),
                       xlab='DustPedia (MJy/sr)', ylab='SINGS (MJy/sr)',
                       figsize=(8,8), title='M82 '+phot)
            p0.save(path_cal+'SINGS-DP_'+phot)
            
            ## IRC (zone 0) - SINGS
            p = pplot(fmt='.', clib='xkcd', legend='upper left',
                      xlog=1, ylog=1, nonposx='clip', nonposy='clip',
                      xlim=(1e-2,1e3), ylim=(1e-2,1e3),
                      xlab='IRC (MJy/sr)', ylab='SINGS (MJy/sr)', 
                      figsize=(8,8), title='M82 '+phot+' inter-calibration')
        else:
            ## IRC (zone 1-5) - SINGS
            p.add_plot(pix_spec, pix_phot2, xerr=pix_spec_unc, yerr=pix_phot2_unc,
                       fmt='.', c=p.clib[z], ec='r')

        ## Linear fit
        popt, pcov = curve_fit(f_lin0, pix_spec[mask], pix_phot2[mask],
                               sigma=pix_phot2_unc[mask])
        calib_fac = popt[0]
        # calib_off = popt[1]
        # popt, pcov = curve_fit(f_lin0, pix_phot2[mask], pix_spec[mask])
        # calib_fac = 1. / popt[0]
        # calib_off = -popt[1] / popt[0]
        
        print('Inter-Calibration ('+phot+') factor ('+zname+') = {:.4}'.format(calib_fac))
        # print('Inter-Calibration ('+phot+') offset ('+zname+') = {:.4}'.format(calib_off))
        label = 'calib fac ('+zname+') = {:.4}'.format(calib_fac)
        # label = zname+': y={0:.4}x+{1:.4}'.format(calib_fac, calib_off)

        p.add_plot(xgrid, f_lin0(xgrid, *popt),
                   c=p.clib[z], ls='-', label=label)

    elif z==Nz-1:
        ## IRC (zone 6) - SINGS (no SINGS data, use zone 5 factor)
        print('Inter-Calibration ('+phot+') factor ('+zname+') = {:.4}'.format(calib_fac))
        # print('Inter-Calibration ('+phot+') offset ('+zname+') = {:.4}'.format(calib_off))
    
    p.set_font()
    
    p.save(path_cal+'calib_'+phot)

    
    ## Spectral correction
    ##---------------------
    if z==0:
        do_correct = input("Do IRC spectral correction (y/n)? ")
        
        sc0 = intercalib(out_irc+'_0')
        sc0.specorrect(filOUT=out_irc)
    else:
        if do_correct:
            sc = intercalib(out_irc)
            sc.specorrect(filOUT=out_irc,
                          factor=calib_fac, offset=0,#calib_off,
                          xlim=(xmin,xmax), ylim=(ymin,ymax))

##---------------------------
##       Plot spectra
##---------------------------
do_plot = input("Plot IRC spectra (y/n)? ")
if do_plot:
    ds = read_fits(out_irc, out_irc+'_unc')
    # ds.data = respect().smooth(out_irc, out_irc, lim_unc=1.e2, cmin=2)
    ds0 = read_fits(out_irc+'_0')
    wcen = intercalib().wcenter(phot)
    for x in range(Nx):
        for y in range(Ny):
            if ~np.isnan(ds.data[:,y,x]).any():
                p = pplot(ds.wave, ds.data[:,y,x], yerr=ds.unc[:,y,x],
                          # xlog=1, ylog=1,
                          c='k', lw=.7, ec='r', legend='upper left')
                p.add_plot(wcen[0], data_phot2[y,x], yerr=unc_phot2[y,x],
                           c='m', marker='o', ms=10, zorder=100, label=phot)
                p.add_plot(wcen[0], data_spec[y,x], yerr=unc_spec[y,x],
                           c='g', marker='^', ms=10, zorder=101, label='IRC-'+phot)
                p.add_plot(ds0.wave, ds0.data[:,y,x],
                           c='y', lw=.7, ls='--', zorder=-1)
                p.save(path_fig+'IRC_('+str(x)+','+str(y)+')')
