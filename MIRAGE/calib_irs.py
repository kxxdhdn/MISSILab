#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""

Inter-calibrate Spitzer/IRS spectra with Spitzer/IRAC4 & MIPS1

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
from laputan.imaging import ( iconvolve, iswarp, iuncert,
                              respect, Jy_per_pix_to_MJy_per_sr )
from laputan.astrom import fixwcs
from laputan.calib import intercalib
from laputan.maths import f_lin, f_lin0
from laputan.plots import pplot

##---------------------------
##       Preliminaries
##---------------------------
## Local
from buildinfo import ( src, Nmc, verbose,
                        path_phot, path_ker, path_cal,
                        path_idl, csv_ker, path_tmp,
                        out_irs, path_fig
)
# Nmc = 2

## Reprojection footprint
refheader = fixwcs(out_irs+'_0'+fitsext).header
swp = iswarp(refheader=refheader,
             tmpdir=path_tmp, verbose=verbose)

Nx = refheader['NAXIS1']
Ny = refheader['NAXIS2']

phot = ['IRAC4', 'MIPS1']
Nphot = len(phot)

## Prepare photometry
##--------------------
precalib = input("Run inter-calibration prepipeline (y/n)? ")

for i in range(Nphot):
    
    if phot[i]=='IRAC4':
        phot_ker = path_ker+'Kernel_HiRes_IRAC_8.0_to_Gauss_06.0'
    elif phot[i]=='MIPS1':
        phot_ker = path_ker+'Kernel_HiRes_MIPS_24_to_Gauss_06.0'
    
    raw_phot1 = path_phot+src+'_'+phot[i]+'_DP'
    fits_phot1 = path_cal+src+'_'+phot[i]+'_DP'
    
    raw_phot2 = path_phot+src+'_'+phot[i]+'_SINGS'
    wgt_phot2 = path_phot+src+'_'+phot[i]+'_SINGS_wgt'
    fits_phot2 = path_cal+src+'_'+phot[i]+'_SINGS'
    
    fits_spec = path_cal+src+'_'+phot[i]+'_IRS'
    
    if precalib:
    
        ## DustPedia (phot1)
        ##===================
        
        if i==0:
            ## Convert phot unit (Jy/pix -> MJy/sr)
            ## This step should be before iswarp reprojection (pixscale changed)
            Jy_per_pix_to_MJy_per_sr(raw_phot1, filOUT=fits_phot1)
            Jy_per_pix_to_MJy_per_sr(raw_phot1+'_unc', filOUT=fits_phot1+'_unc')

            for j in trange(Nmc+1,# leave=False,
                    desc='Generating DustPedia uncertainty map [MC]'):
                if j==0:
                    ## Reproject phot
                    comb = swp.combine(fits_phot1,# combtype='wgt_avg',
                                       cropedge=False, tmpdir=path_tmp,
                                       filOUT=fits_phot1).image
                    comb = swp.combine(fits_phot1+'_unc',# combtype='wgt_avg',
                                       cropedge=False, tmpdir=path_tmp,
                                       filOUT=fits_phot1+'_unc').image

                    ## Convolve phot
                    conv = iconvolve(fits_phot1, phot_ker, csv_ker, filOUT=fits_phot1)
                else:
                    comb = swp.combine(fits_phot1,# combtype='wgt_avg',
                                       cropedge=False, tmpdir=path_tmp, uncpdf='norm',
                                       filOUT=path_tmp+src+'_'+phot[i]+'_DP_'+str(j)).image
                    
                    conv = iconvolve(path_tmp+src+'_'+phot[i]+'_DP_'+str(j),
                                     kfile=phot_ker, klist=csv_ker,
                                     filOUT=path_tmp+src+'_'+phot[i]+'_DP_'+str(j))
            
                conv.do_conv(idldir=path_idl)
            
            ## Propagate DustPedia IRAC4 uncertainty
            mcimage = []
            for j in range(Nmc+1):
                if j==0:
                    hdr = read_fits(fits_phot1).header
                else:
                    mcimage.append(read_fits(path_tmp+src+'_'+phot[i]+'_DP_'+str(j)).data)
            if Nmc>1:
                mcimage = np.array(mcimage)
                unc = np.nanstd(mcimage, axis=0)
                write_fits(fits_phot1+'_unc', hdr, unc)
        
        ## SINGS (phot2)
        ##===============
        
        ## Create uncertainty via weight map (suppose uniform contribution)
        ## The weight maps contain the information on the number of frames 
        ## that were used to create the science mosaics at each pixel (value= # frames x10); 
        ## the pixel size of the weight maps is the same as the science mosaics.
        ## -- SINGS v5 release doc
        if phot[i]=='IRAC4':
            bg = read_fits(raw_phot2).data[502:542,1067:1107]
            bg_wgt = read_fits(wgt_phot2).data[502:542,1067:1107] / 10.
        elif phot[i]=='MIPS1':
            bg = read_fits(raw_phot2).data[450:490,1606:1646]
            bg_wgt = read_fits(wgt_phot2).data[450:490,1606:1646] / 10.
        iuncert(raw_phot2, filOUT=raw_phot2+'_unc',
                filWGT=wgt_phot2, wfac=.1,
                BG_image=bg, BG_weight=bg_wgt)
        
        for j in trange(Nmc+1,# leave=False,
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
                                   filOUT=path_tmp+src+'_'+phot[i]+'_SINGS_'+str(j)).image
            
                conv = iconvolve(path_tmp+src+'_'+phot[i]+'_SINGS_'+str(j),
                                 kfile=phot_ker, klist=csv_ker,
                                 filOUT=path_tmp+src+'_'+phot[i]+'_SINGS_'+str(j))
            conv.do_conv(idldir=path_idl)
        
        mcimage = []
        for j in range(Nmc+1):
            if j==0:
                hdr = read_fits(fits_phot2).header
            else:
                mcimage.append(read_fits(path_tmp+src+'_'+phot[i]+'_SINGS_'+str(j)).data)
        if Nmc>1:
            mcimage = np.array(mcimage)
            unc = np.nanstd(mcimage, axis=0)
            write_fits(fits_phot2+'_unc', hdr, unc)

        ## Synthetic photometry (spec)
        ##=============================
        mcimage = []
        for j in trange(Nmc+1,# leave=False,
                        desc='IRS sythetic photometry [MC]'):
            ic = intercalib(out_irs+'_'+str(j))
            sp = ic.synthetic_photometry(phot[i])
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
    if i==0:
        data_phot1 = read_fits(fits_phot1).data
        unc_phot1 = read_fits(fits_phot1+'_unc').data
        pix_phot1 = data_phot1[:,:].reshape((-1,))
        
    data_phot2 = read_fits(fits_phot2).data
    unc_phot2 = read_fits(fits_phot2+'_unc').data
    data_spec = read_fits(fits_spec).data
    unc_spec = read_fits(fits_spec+'_unc').data
    if i==0:
        d1_phot2 = data_phot2
        u1_phot2 = unc_phot2
        d1_spec = data_spec
        u1_spec = unc_spec
    elif i==1:
        d2_phot2 = data_phot2
        u2_phot2 = unc_phot2
        d2_spec = data_spec
        u2_spec = unc_spec
    
    ## Plot
    ##------
    pix_phot2 = data_phot2[:,:].reshape((-1,))
    pix_spec = data_spec[:,:].reshape((-1,))
    pix_phot2_unc = unc_phot2[:,:].reshape((-1,))
    pix_spec_unc = unc_spec[:,:].reshape((-1,))
    
    # xgrid = np.arange(0,1e3,1)
    xgrid = np.logspace(-3,3,10000)
    
    ## NaNs mask
    mask = ~np.ma.array(pix_spec,
                        mask=np.logical_or(
                            np.isnan(pix_spec),np.isnan(pix_phot2))).mask

    # if i==0:
    #     mask_all = mask
    # else:
    #     mask_all = np.ma.array(mask,
    #                            mask=np.logical_and(mask,mask_all)).mask
    
    ## S/N ratio
    # print('S/N ('+phot[i]+' - IRS) = \n', pix_spec[mask]/pix_spec_unc[mask])
    # print('S/N ('+phot[i]+' - SINGS) = \n', pix_phot2[mask]/pix_phot2_unc[mask])
    # exit()
        
    if i==0:
        ## SINGS - DP
        p0 = pplot(pix_phot1,  pix_phot2, yerr=pix_phot2_unc,
                   fmt='.', c='y', ec='r',elw=10,
                   xlog=1, ylog=1,
                   xlim=(0,1e3), ylim=(0,1e3),
                   xlab='DustPedia (MJy/sr)', ylab='SINGS (MJy/sr)',
                   figsize=(8,8), title='M82 IRAC4')
        p0.save(path_cal+'SINGS-DP_IRAC4')
        
    ## IRS - SINGS
    p = pplot(fmt='.', legend='upper left',
              xlog=1, ylog=1,
              xlim=(0,1e3), ylim=(0,1e3),
              xlab='IRS (MJy/sr)', ylab='SINGS (MJy/sr)', 
              figsize=(8,8), title='M82 '+phot[i]+' inter-calibration')
    ## IRS - SINGS
    p.add_plot(pix_spec, pix_phot2, xerr=pix_spec_unc, yerr=pix_phot2_unc,
               fmt='.', c='y', ec='r')
    
    ## Linear fit
    popt, pcov = curve_fit(f_lin, pix_spec[mask], pix_phot2[mask],
                           sigma=pix_spec_unc[mask])
    calib_fac = popt[0]
    calib_off = popt[1]
    # popt, pcov = curve_fit(f_lin0, pix_phot2[mask], pix_spec[mask])
    # calib_fac = 1. / popt[0]
    # calib_off = -popt[1] / popt[0]
    
    print('Inter-Calibration ('+phot[i]+') factor = {:.4}'.format(calib_fac))
    print('Inter-Calibration ('+phot[i]+') offset = {:.4}'.format(calib_off))
    # label = phot[i]+' calib fac = {:.4}'.format(calib_fac)
    label = src+': y={0:.4}x+{1:.4}'.format(calib_fac, calib_off)
    
    p.add_plot(xgrid, f_lin(xgrid, *popt),
               c='k', ls='-', label=label)
    
    p.set_font()
    
    p.save(path_cal+'calib_'+phot[i])
    
    
    ## Spectral correction
    ##---------------------
    do_correct = input("Do IRS spectral correction [{}] (y/n)? ".format(phot[i]))
    if i==0:
        sc0 = intercalib(out_irs+'_0')
        sc0.specorrect(filOUT=out_irs)
        wmin = 5.21
        wmax = 14.28
    elif i==1:
        wmin = 14.29
        wmax = 38.00
        
    if do_correct:
        sc = intercalib(out_irs)
        sc.specorrect(filOUT=out_irs,
                      factor=calib_fac, offset=calib_off,
                      wlim=(wmin,wmax))

##---------------------------
##       Plot spectra
##---------------------------
do_plot = input("Plot IRS spectra (y/n)? ")
if do_plot:
    ds = read_fits(out_irs, out_irs+'_unc')
    # ds.data = respect().smooth(out_irs, out_irs, lim_unc=1.e2, cmin=2)
    ds0 = read_fits(out_irs+'_0')
    wcen = intercalib().wcenter(phot)
    for x in range(Nx):
        for y in range(Ny):
            if ~np.isnan(ds.data[:,y,x]).any():
                p = pplot(ds.wave, ds.data[:,y,x], yerr=ds.unc[:,y,x],
                          c='k', lw=.7, ec='r', xlog=1, ylog=1, legend='upper left')
                p.add_plot(wcen[0], d1_phot2[y,x], yerr=u1_phot2[y,x],
                           c='m', marker='o', ms=10, zorder=100, label=phot[0])
                p.add_plot(wcen[0], d1_spec[y,x], yerr=u1_spec[y,x],
                           c='g', marker='^', ms=10, zorder=101, label='IRS-'+phot[0])
                p.add_plot(wcen[1], d2_phot2[y,x], yerr=u1_phot2[y,x],
                           c='orange', marker='o', ms=10, zorder=102, label=phot[1])
                p.add_plot(wcen[1], d2_spec[y,x], yerr=u1_spec[y,x],
                           c='c', marker='^', ms=10, zorder=103, label='IRS-'+phot[1])
                p.add_plot(ds0.wave, ds0.data[:,y,x],
                           c='y', ls='--', zorder=-1)
                p.save(path_fig+'IRS_('+str(x)+','+str(y)+')')
