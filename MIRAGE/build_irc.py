#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""

Build AKARI/IRC spectral cube
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
from laputan.imaging import ( iconvolve, iswarp, sextract, irebin, iuncert,
                              imontage, improve, Jy_per_pix_to_MJy_per_sr )
from laputan.astrom import fixwcs
from laputan.calib import intercalib
from laputan.maths import f_lin, f_lin0
from laputan.plots import pplot


##---------------------------
##       Preliminaries
##---------------------------
## Local
from buildinfo import ( src, Nmc, verbose, coadd_tool,
                        path_irc, path_build, parobs, # build
                        fits_irc, out_irc,
                        path_phot, path_ker, path_cal, # calib
                        path_idl, csv_ker, path_tmp, path_fig
)
Nobs = len(parobs)
Nmc = 20

##---------------------------
##        Build slits
##---------------------------
do_build = input("(Re)build IRC slits (y/n)? ")
if do_build=='y':

    slit_length = 24
    slit_width = 1
    Nsub = 6
    for i in trange(Nobs, #leave=False,
                    desc='<sextract> IRC slit building'):
    
        ## Reproduce Y12 spectra
        ##-----------------------
        # if fits_irc[i][-1]=='s':
        #     Nsub = 2 # Ns
        # if fits_irc[i][-1]=='h':
        #     Nsub = 6 # Nh
        
        ## MC add unc
        for j in trange(Nmc+1, leave=False,
                        desc='IRC slit unc [MC]'):
            sext = sextract(path_irc, parobs[i])
            if j==0:
                sext.spec_build(fits_irc[i],
                                Nx=slit_width, Ny=slit_length, Nsub=Nsub)
            else:
                imp = improve(fits_irc[i])
                imp.rand_splitnorm([fits_irc[i]+'_unc_N',fits_irc[i]+'_unc_P'])
                write_fits(fits_irc[i]+'_'+str(j),
                           imp.hdr, imp.im, imp.wvl, imp.wmod)

##---------------------------
##    Prepare photometry
##---------------------------
phot = 'IRAC1'
phot_ker = path_ker+'Kernel_HiRes_IRAC_3.6_to_Gauss_06.0'

if not os.path.exists(path_tmp+'calib/'):
    os.makedirs(path_tmp+'calib/')
    
raw_phot1 = path_phot+src+'_'+phot+'_DP'
fits_phot1 = path_cal+src+'_'+phot+'_DP'
tmp_phot1 = path_tmp+'calib/'+src+'_'+phot+'_DP'

raw_phot2 = path_phot+src+'_'+phot+'_SINGS'
wgt_phot2 = path_phot+src+'_'+phot+'_SINGS_wgt'
fits_phot2 = path_cal+src+'_'+phot+'_SINGS'
tmp_phot2 = path_tmp+'calib/'+src+'_'+phot+'_SINGS'

fits_spec = path_cal+src+'_'+phot+'_IRC'
tmp_spec = path_tmp+'calib/'+src+'_'+phot+'_IRC'

pre_calib1 = input("Run inter-calibration prepipeline (DustPedia) (y/n)? ")
pre_calib2 = input("Run inter-calibration prepipeline (SINGS) (y/n)? ")

## DustPedia (phot1)
##===================
if pre_calib1=='y':
    
    ## Convert phot unit (Jy/pix -> MJy/sr)
    ## This step should be before reprojection (pixscale may change)
    Jy_per_pix_to_MJy_per_sr(raw_phot1, filOUT=tmp_phot1)

    ## DustPedia IRAC1 of M82 has no uncertainty map
    ## Create homogeneous uncertainty map
    bg = read_fits(tmp_phot1).data[1713:1753,3110:3150]
    iuncert(tmp_phot1, filOUT=tmp_phot1+'_unc', BG_image=bg)
    
    for f in fits_irc:
        ## Reproject phot
        fname = os.path.basename(f)
        refheader = fixwcs(f+fitsext).header
        mtg = imontage('exact', tmpdir=path_tmp)
        mtg.reproject_mc(tmp_phot1, refheader=refheader,
                         dist='norm', Nmc=Nmc,
                         filOUT=tmp_phot1+'_'+fname)

        ## Convolve phot
        for j in trange(Nmc+1, leave=False,
                        desc=fname+': DP conv [MC]'):
            if j==0:
                conv = iconvolve(tmp_phot1+'_'+fname,
                                 kfile=phot_ker, klist=csv_ker,
                                 filOUT=fits_phot1+'_'+fname)
            else:
                conv = iconvolve(tmp_phot1+'_'+fname+'_'+str(j),
                                 kfile=phot_ker, klist=csv_ker,
                                 # filUNC=tmp_phot1+'_unc', dist='norm',
                                 filOUT=tmp_phot1+'_'+fname+'_'+str(j))
            conv.do_conv(idldir=path_idl)

        ## Calculate unc
        mcimage = []
        for j in range(Nmc+1):
            if j==0:
                hdr = read_fits(fits_phot1+'_'+fname).header
            else:
                mcimage.append(read_fits(tmp_phot1+'_'+fname+'_'+str(j)).data)
        if Nmc>1:
            mcimage = np.array(mcimage)
            unc = np.nanstd(mcimage, axis=0)
            write_fits(fits_phot1+'_'+fname+'_unc', hdr, unc)

## SINGS (phot2)
##===============
if pre_calib2=='y':
    
    ## Create uncertainty via weight map (suppose uniform contribution)
    ## The weight maps contain the information on the number of frames 
    ## that were used to create the science mosaics at each pixel (value= # frames x10); 
    ## the pixel size of the weight maps is the same as the science mosaics.
    ## -- SINGS v5 release doc
    # gen_unc2 = input("Create SINGS IRAC1 uncertainty map (y/n)? ")
    # if gen_unc2=='y':
    #     bg = read_fits(raw_phot2).data[694:734,56:96]
    #     bg_wgt = read_fits(wgt_phot2).data[694:734,56:96] / 10.
    #     iuncert(raw_phot2, filOUT=raw_phot2+'_unc',
    #             filWGT=wgt_phot2, wfac=.1,
    #             BG_image=bg, BG_weight=bg_wgt)

    for f in fits_irc:
        ## Reproject phot
        fname = os.path.basename(f)
        refheader = fixwcs(f+fitsext).header
        mtg = imontage('exact', tmpdir=path_tmp)
        mtg.reproject_mc(raw_phot2, refheader=refheader,
                         dist='norm', Nmc=Nmc,
                         filOUT=tmp_phot2+'_'+fname)

        ## Convolve phot
        for j in trange(Nmc+1, leave=False,
                        desc=fname+': SINGS conv [MC]'):
            if j==0:
                conv = iconvolve(tmp_phot2+'_'+fname,
                                 kfile=phot_ker, klist=csv_ker,
                                 filOUT=fits_phot2+'_'+fname)
            else:
                conv = iconvolve(tmp_phot2+'_'+fname+'_'+str(j),
                                 kfile=phot_ker, klist=csv_ker,
                                 # filUNC=tmp_phot2+'_unc', dist='norm',
                                 filOUT=tmp_phot2+'_'+fname+'_'+str(j))
            conv.do_conv(idldir=path_idl)
        
        ## Calculate unc
        mcimage = []
        for j in range(Nmc+1):
            if j==0:
                hdr = read_fits(fits_phot2+'_'+fname).header
            else:
                mcimage.append(read_fits(tmp_phot2+'_'+fname+'_'+str(j)).data)
        if Nmc>1:
            mcimage = np.array(mcimage)
            unc = np.nanstd(mcimage, axis=0)
            write_fits(fits_phot2+'_'+fname+'_unc', hdr, unc)
        
        ## Synthetic photometry (spec)
        ##=============================
        mcimage = []
        for j in trange(Nmc+1, leave=False,
                        desc=fname+': IRC sythetic photometry [MC]'):
            if j==0:
                ic = intercalib(f)
                sp = ic.synthetic_photometry(phot)
                write_fits(fits_spec+'_'+fname, ic.hdr, sp.Fnu_filt)
            else:
                ic = intercalib(f+'_'+str(j))
                sp = ic.synthetic_photometry(phot)
                ## Cal unc
                mcimage.append(sp.Fnu_filt)
        if Nmc>1:
            # print('Calculating uncertainty cube...')
            mcimage = np.array(mcimage)
            unc = np.nanstd(mcimage, axis=0)
            write_fits(fits_spec+'_'+fname+'_unc', ic.hdr, unc)

##---------------------------
##     Inter-Calibration
##---------------------------
int_calib = input("Run inter-calibration (y/n)? ")
if int_calib=='y':
    do_correct = input("Do IRC spectral correction (y/n)? ")
    
    for f in fits_irc:
        fname = os.path.basename(f)

        ## Data
        ##------
        data_phot1 = read_fits(fits_phot1+'_'+fname).data
        unc_phot1 = read_fits(fits_phot1+'_'+fname+'_unc').data
        data_phot2 = read_fits(fits_phot2+'_'+fname).data
        unc_phot2 = read_fits(fits_phot2+'_'+fname+'_unc').data
        data_spec = read_fits(fits_spec+'_'+fname).data
        unc_spec = read_fits(fits_spec+'_'+fname+'_unc').data
            
        ## Plot
        ##------
        pix_phot1 = data_phot1[:,:].reshape((-1,))
        pix_phot1_unc = unc_phot1[:,:].reshape((-1,))
        pix_phot2 = data_phot2[:,:].reshape((-1,))
        pix_phot2_unc = unc_phot2[:,:].reshape((-1,))
        pix_spec = data_spec[:,:].reshape((-1,))
        pix_spec_unc = unc_spec[:,:].reshape((-1,))
        
        # xgrid = np.arange(0,1e3,1)
        xgrid = np.logspace(-3,3,10000)
        
        ## NaN mask
        mask = ~np.ma.array(pix_spec,
                            mask=np.logical_or(
                                np.isnan(pix_spec),np.isnan(pix_phot2))).mask
        
        ## S/N ratio
        # print('S/N (IRC) = \n', pix_spec[mask]/pix_spec_unc[mask])
        # print('S/N (SINGS) = \n', pix_phot2[mask]/pix_phot2_unc[mask])
        # exit()

        if not '5124077' in fname:
            calib_fac = data_phot2/data_spec
            
            ## Spectral correction
            ##---------------------            
            if do_correct=='y':
                for j in trange(Nmc+1, leave=False,
                                desc=fname+': spectral correction [MC]'):
                    if j==0:
                        imp = improve(f)
                    else:
                        imp = improve(f+'_'+str(j))
                    im = imp.im * calib_fac
                    write_fits(f+'_'+str(j), imp.hdr, im, imp.wvl, imp.wmod)
            else:
                sc = intercalib(f)
                sc.specorrect(filOUT=f+'_0')
        else:
            sc = intercalib(f)
            sc.specorrect(filOUT=f+'_0')
        
        '''
        if not '5124077' in fname:
            ## SINGS - DP
            p0 = pplot(pix_phot1,  pix_phot2,
                       yerr=pix_phot2_unc, xerr=pix_phot1_unc,
                       fmt='.', c='k', ec='r', elw=1,
                       xlog=1, ylog=1,
                       # xlim=(0,1e3), ylim=(0,1e3),
                       xlab='DustPedia (MJy/sr)', ylab='SINGS (MJy/sr)',
                       figsize=(8,8), title=src+'_'+phot+'_'+fname)
            p0.save(path_cal+'SINGS-DP_'+phot+'_'+fname+'.png')
            
            ## IRC - SINGS
            p = pplot(pix_phot2, pix_spec,
                      yerr=pix_spec_unc, xerr=pix_phot2_unc,
                      fmt='.', c='k', ec='r', elw=1, legend='upper left',
                      xlog=1, ylog=1, nonposx='clip', nonposy='clip',
                      # xlim=(1e-2,1e3), ylim=(1e-2,1e3),
                      xlab='SINGS (MJy/sr)', ylab='IRC (MJy/sr)', 
                      figsize=(8,8), title=src+'_'+phot+'_'+fname)

            ## Linear fit
            # popt, pcov = curve_fit(f_lin, pix_spec[mask], pix_phot2[mask],
            #                        sigma=pix_phot2_unc[mask])
            # calib_fac = popt[0]
            # calib_off = popt[1]
            
            popt, pcov = curve_fit(f_lin, pix_phot2[mask], pix_spec[mask],
                                   sigma=pix_spec_unc[mask])
            calib_fac = 1. / popt[0]
            calib_off = -popt[1] / popt[0]
            
            print('Inter-Calibration ('+phot+') factor ('+fname+') = {:.4}'.format(calib_fac))
            print('Inter-Calibration ('+phot+') offset ('+fname+') = {:.4}'.format(calib_off))
            # label = 'calib fac ('+fname+') = {:.4}'.format(calib_fac)
            label = fname+': y={0:.4}x+{1:.4}'.format(calib_fac, calib_off)
            
            p.add_plot(xgrid, f_lin(xgrid, *popt),
                       c='y', ls='-', label=label)
            
            p.set_font()
            
            p.save(path_cal+'intercalib_'+phot+'_'+fname+'.png')
            
            ## Spectral correction
            ##---------------------            
            if do_correct:
                for j in trange(Nmc+1, leave=False,
                                desc=fname+': spectral correction [MC]'):
                    if j==0:
                        sc = intercalib(f)
                        sc.specorrect(filOUT=f+'_0', factor=calib_fac, offset=calib_off)
                    else:
                        sc = intercalib(f+'_'+str(j))
                        sc.specorrect(filOUT=f+'_'+str(j), factor=calib_fac, offset=calib_off)
            else:
                sc = intercalib(f)
                sc.specorrect(filOUT=f+'_0')
        else:
            sc = intercalib(f)
            sc.specorrect(filOUT=f+'_0')
        '''
##---------------------------
##        Coadd slits
##---------------------------

## <iswarp> coadding
##===================
swp_coadd = input("Coadd IRC slits (y/n)? ")
# swp_coadd = None
if swp_coadd=='y':
    refheader = fixwcs(fits_irc[0]+fitsext).header
    swp = iswarp(fits_irc, refheader=refheader, 
                 # center='9:55:52,69:40:45', pixscale=6,
                 tmpdir=path_build, verbose=verbose)
    
    ## Reprendre MC
    ##--------------
    for j in trange(Nmc+1, #leave=False, 
                    desc='<iswarp> IRC Coadding [MC]'):

        ## Uncalibrated
        # comb = swp.combine(fits_irc, combtype='wgt_avg',
        #                    keepedge=True, cropedge=True,
        #                    tmpdir=path_build+'MC_'+str(j)+'/',
        #                    filOUT=out_irc)
        # irebin(out_irc, filOUT=out_irc, pixscale=6,
        #        total=False, extrapol=True, verbose=True)

        ## Calibrated
        fits_irc_mc = []
        for f in fits_irc:
            fits_irc_mc.append(f+'_'+str(j))
        comb = swp.combine(fits_irc_mc, #combtype='wgt_avg',
                           keepedge=True, cropedge=True,
                           tmpdir=path_build+'MC_'+str(j)+'/',
                           filOUT=out_irc+'_'+str(j))
        ## Rescale pixel size
        irebin(out_irc+'_'+str(j), filOUT=out_irc+'_'+str(j), pixscale=6,
               total=False, extrapol=True, verbose=True)
        
    ## Cal unc
    ##---------
    mcimage = []
    for j in trange(Nmc+1, leave=False,
                    desc='IRC Reading [MC]'):
        hd = read_fits(out_irc+'_'+str(j))
        if j==0:
            header = hd.header
            wvl = hd.wave
        else:
            mcimage.append(hd.data)
    if Nmc>1:
        print('Calculating uncertainty cube...')
        mcimage = np.array(mcimage)
        unc = np.nanstd(mcimage, axis=0)
        write_fits(out_irc+'_unc', header, unc, wvl)

## <imontage> coadding
##=====================
# mtg_coadd = input("Coadd IRC slits (y/n)? ")
mtg_coadd = None
if mtg_coadd=='y':
    refheader = fixwcs(fits_irc[0]+fitsext).header # one slit frame
    swp = iswarp(fits_irc, refheader=refheader, 
                 tmpdir=path_build, verbose=verbose)
    refheader = swp.refheader # expanded frame
    
    mtg = imontage('exact', tmpdir=path_tmp)

    ## Uncalibrated
    # mtg.coadd(fits_irc, refheader=refheader,
    #           filOUT=out_irc)
    # irebin(out_irc, filOUT=out_irc, pixscale=6,
    #        total=False, extrapol=True, verbose=True)

    ## Calibrated
    for j in trange(Nmc+1, leave=False,
                    desc='<imontage> IRC coadding [MC]'):
        fits_irc_mc = []
        for f in fits_irc:
            fits_irc_mc.append(f+'_'+str(j))
        mtg.coadd(fits_irc_mc, refheader=refheader,
                  filOUT=out_irc+'_'+str(j))
        ## Rescale pixel size
        irebin(out_irc+'_'+str(j), filOUT=out_irc+'_'+str(j), pixscale=6,
               total=False, extrapol=True, verbose=True)
        
    ## Cal unc
    ##---------
    mcimage = []
    for j in trange(Nmc+1, leave=False,
                    desc='IRC Reading [MC]'):
        hd = read_fits(out_irc+'_'+str(j))
        if j==0:
            header = hd.header
            wvl = hd.wave
        else:
            mcimage.append(hd.data)
    if Nmc>1:
        print('Calculating uncertainty cube...')
        mcimage = np.array(mcimage)
        unc = np.nanstd(mcimage, axis=0)
        write_fits(out_irc+'_unc', header, unc, wvl)
        
##---------------------------
##  Atlas inter-Calibration
##---------------------------
pre_calib1 = input("Run atlas inter-calibration prepipeline (DustPedia) (y/n)? ")
pre_calib2 = input("Run atlas inter-calibration prepipeline (SINGS) (y/n)? ")

## DustPedia (phot1)
##===================
if pre_calib1=='y':
    
    ## Convert phot unit (Jy/pix -> MJy/sr)
    ## This step should be before reprojection (pixscale may change)
    Jy_per_pix_to_MJy_per_sr(raw_phot1, filOUT=path_tmp+src+'_'+phot+'_DP')

    ## DustPedia IRAC1 of M82 has no uncertainty map
    ## Create homogeneous uncertainty map
    bg = read_fits(path_tmp+src+'_'+phot+'_DP').data[1713:1753,3110:3150]
    iuncert(path_tmp+src+'_'+phot+'_DP',
            filOUT=path_tmp+src+'_'+phot+'_DP_unc', BG_image=bg)
    
    ## Reproject phot
    refheader = fixwcs(out_irc+'_0'+fitsext).header
    mtg = imontage('exact', tmpdir=path_tmp)
    mtg.reproject_mc(path_tmp+src+'_'+phot+'_DP', refheader=refheader,
                     dist='norm', Nmc=Nmc, filOUT=tmp_phot1)

    ## Convolve phot
    for j in trange(Nmc+1, leave=False,
                    desc=phot+': DP conv [MC]'):
        if j==0:
            conv = iconvolve(tmp_phot1,
                             kfile=phot_ker, klist=csv_ker,
                             filOUT=fits_phot1)
        else:
            conv = iconvolve(tmp_phot1+'_'+str(j),
                             kfile=phot_ker, klist=csv_ker,
                             # filUNC=fits_phot1+'_unc', dist='norm',
                             filOUT=tmp_phot1+'_'+str(j))
        conv.do_conv(idldir=path_idl)

    ## Calculate unc
    mcimage = []
    for j in range(Nmc+1):
        if j==0:
            hdr = read_fits(fits_phot1).header
        else:
            mcimage.append(read_fits(tmp_phot1+'_'+str(j)).data)
    if Nmc>1:
        mcimage = np.array(mcimage)
        unc = np.nanstd(mcimage, axis=0)
        write_fits(fits_phot1+'_unc', hdr, unc)

## SINGS (phot2)
##===============
if pre_calib2=='y':
    
    ## Create uncertainty via weight map (suppose uniform contribution)
    ## The weight maps contain the information on the number of frames 
    ## that were used to create the science mosaics at each pixel (value= # frames x10); 
    ## the pixel size of the weight maps is the same as the science mosaics.
    ## -- SINGS v5 release doc
    # gen_unc2 = input("Create SINGS IRAC1 uncertainty map (y/n)? ")
    # if gen_unc2=='y':
    #     bg = read_fits(raw_phot2).data[694:734,56:96]
    #     bg_wgt = read_fits(wgt_phot2).data[694:734,56:96] / 10.
    #     iuncert(raw_phot2, filOUT=raw_phot2+'_unc',
    #             filWGT=wgt_phot2, wfac=.1,
    #             BG_image=bg, BG_weight=bg_wgt)

    ## Reproject phot
    refheader = fixwcs(out_irc+'_0'+fitsext).header
    mtg = imontage('exact', tmpdir=path_tmp)
    mtg.reproject_mc(raw_phot2, refheader=refheader,
                     dist='norm', Nmc=Nmc,
                     filOUT=tmp_phot2)

    ## Convolve phot
    for j in trange(Nmc+1, leave=False,
                    desc=phot+': SINGS conv [MC]'):
        if j==0:
            conv = iconvolve(tmp_phot2,
                             kfile=phot_ker, klist=csv_ker,
                             filOUT=fits_phot2)
        else:
            conv = iconvolve(tmp_phot2+'_'+str(j),
                             kfile=phot_ker, klist=csv_ker,
                             # filUNC=fits_phot2+'_unc', dist='norm',
                             filOUT=tmp_phot2+'_'+str(j))
        conv.do_conv(idldir=path_idl)
    
    ## Calculate unc
    mcimage = []
    for j in range(Nmc+1):
        if j==0:
            hdr = read_fits(fits_phot2).header
        else:
            mcimage.append(read_fits(tmp_phot2+'_'+str(j)).data)
    if Nmc>1:
        mcimage = np.array(mcimage)
        unc = np.nanstd(mcimage, axis=0)
        write_fits(fits_phot2+'_unc', hdr, unc)
    
    ## Synthetic photometry (spec)
    ##=============================
    mcimage = []
    for j in trange(Nmc+1, leave=False,
                    desc=phot+': IRC sythetic photometry [MC]'):
        ic = intercalib(out_irc+'_'+str(j))
        sp = ic.synthetic_photometry(phot)
        if j==0:
            write_fits(fits_spec, ic.hdr, sp.Fnu_filt)
        else:
            ## Cal unc
            mcimage.append(sp.Fnu_filt)
    if Nmc>1:
        # print('Calculating uncertainty cube...')
        mcimage = np.array(mcimage)
        unc = np.nanstd(mcimage, axis=0)
        write_fits(fits_spec+'_unc', ic.hdr, unc)

## Data
##------
data_phot1 = read_fits(fits_phot1).data
unc_phot1 = read_fits(fits_phot1+'_unc').data
data_phot2 = read_fits(fits_phot2).data
unc_phot2 = read_fits(fits_phot2+'_unc').data
data_spec = read_fits(fits_spec).data
unc_spec = read_fits(fits_spec+'_unc').data

all_calib = input("Run atlas inter-calibration (y/n)? ")
if all_calib=='y':
    
    ## Plot
    ##------
    pix_phot1 = data_phot1[:,:].reshape((-1,))
    pix_phot1_unc = unc_phot1[:,:].reshape((-1,))
    pix_phot2 = data_phot2[:,:].reshape((-1,))
    pix_phot2_unc = unc_phot2[:,:].reshape((-1,))
    pix_spec = data_spec[:,:].reshape((-1,))
    pix_spec_unc = unc_spec[:,:].reshape((-1,))
    
    # xgrid = np.arange(0,1e3,1)
    xgrid = np.logspace(-3,3,10000)
    
    ## NaNs mask
    mask = ~np.ma.array(pix_spec,
                        mask=np.logical_or(
                            np.isnan(pix_spec),np.isnan(pix_phot2))).mask
    ## SINGS - DP
    p0 = pplot(pix_phot1,  pix_phot2,
               yerr=pix_phot2_unc, xerr=pix_phot1_unc,
               fmt='.', c='k', ec='r', elw=1,
               xlog=1, ylog=1,
               xlim=(0,1e3), ylim=(0,1e3),
               xlab='DustPedia (MJy/sr)', ylab='SINGS (MJy/sr)',
               figsize=(8,8), title=src+'_'+phot)
    p0.save(path_cal+'SINGS-DP_'+phot+'.png')
            
    ## IRC - SINGS
    p = pplot(pix_phot2, pix_spec,
              yerr=pix_spec_unc, xerr=pix_phot2_unc,
              fmt='.', c='k', ec='r', legend='upper left',
              xlog=1, ylog=1, nonposx='clip', nonposy='clip',
              # xlim=(1e-2,1e3), ylim=(1e-2,1e3),
              xlab='SINGS (MJy/sr)', ylab='IRC (MJy/sr)', 
              figsize=(8,8), title=src+'_'+phot)

    ## Linear fit
    # popt, pcov = curve_fit(f_lin0, pix_spec[mask], pix_phot2[mask],
    #                        sigma=pix_phot2_unc[mask])
    # calib_fac = popt[0]
    # calib_off = popt[1]
    
    popt, pcov = curve_fit(f_lin0, pix_phot2[mask], pix_spec[mask],
                           sigma=pix_spec_unc[mask])
    calib_fac = 1. / popt[0]
    # calib_off = -popt[1] / popt[0]
    
    print('Inter-Calibration ('+phot+') factor = {:.4}'.format(calib_fac))
    # print('Inter-Calibration ('+phot+') offset = {:.4}'.format(calib_off))
    label = 'calib fac = {:.4}'.format(calib_fac)
    # label = phot+': y={0:.4}x+{1:.4}'.format(calib_fac, calib_off)
    
    p.add_plot(xgrid, f_lin0(xgrid, *popt),
               c='y', ls='-', label=label)
    
    p.set_font()
    
    p.save(path_cal+'intercalib_'+phot+'.png')
            
    ## Spectral correction
    ##---------------------    
    sc0 = intercalib(out_irc+'_0')
    sc0.specorrect(filOUT=out_irc)
    
    do_correct = input("Do IRC spectral correction (y/n)? ")
    if do_correct=='y':
        sc = intercalib(out_irc)
        sc.specorrect(filOUT=out_irc, factor=calib_fac, offset=0)
        
##---------------------------
##       Plot spectra
##---------------------------
do_plot = input("Plot IRC spectra (y/n)? ")
if do_plot=='y':
    ds = read_fits(out_irc, out_irc+'_unc')
    # ds.data = respect().smooth(out_irc+'_0', out_irc+'_0', lim_unc=1.e2, cmin=2)
    ds0 = read_fits(out_irc+'_0')
    wcen = intercalib().wcenter(phot)
    Nw, Ny, Nx = ds.data.shape
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

