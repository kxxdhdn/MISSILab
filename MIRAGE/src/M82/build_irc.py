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
## rapyuta
from rapyuta.inout import fclean, fitsext, read_fits, write_fits
from rapyuta.imaging import ( iconvolve, iswarp, sextract, irebin, iuncert,
                              imontage, improve, Jy_per_pix_to_MJy_per_sr )
from rapyuta.astrom import fixwcs
from rapyuta.calib import intercalib
from rapyuta.maths import f_lin, f_lin0, f_lin1
from rapyuta.plots import pplot


##----------------------------------------------------------

##                     Preliminaries

##----------------------------------------------------------

## Local
from buildinfo import ( src, Nmc, verbose, coadd_tool, pixscale,
                        path_irc, path_build, parobs, # build
                        fits_irc, out_irc,
                        path_phot, path_ker, path_cal, # calib
                        path_idl, csv_ker, path_tmp, path_fig
)

## Banner
print('\n============================================================\n')

print('        MIRAGE - AKARI/IRC cube builder - '+src)

print('\n============================================================\n')

Nobs = len(parobs)

# Nmc = 2

##----------------------------------------------------------

##                      Build slits

##----------------------------------------------------------
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

##----------------------------------------------------------

##                   Inter-Calibration

##----------------------------------------------------------
pre_calib1 = input("Run inter-calibration prepipeline (DustPedia) (y/n)? ")
pre_calib2 = input("Run inter-calibration prepipeline (SINGS) (y/n)? ")
int_calib = input("Run inter-calibration (y/n)? ")
do_correct = input("Do IRC spectral correction for each slit (y/n)? ")

if not os.path.exists(path_tmp+'calib/'):
    os.makedirs(path_tmp+'calib/')
    
## Prepare photometry
##--------------------
phot = 'IRAC1'
phot_ker = path_ker+'Kernel_HiRes_IRAC_3.6_to_Gauss_06.0'

raw_phot1 = path_phot+src+'_'+phot+'_DP'
fits_phot1 = path_cal+src+'_'+phot+'_DP'
tmp_phot1 = path_tmp+'calib/'+src+'_'+phot+'_DP'

raw_phot2 = path_phot+src+'_'+phot+'_SINGS'
wgt_phot2 = path_phot+src+'_'+phot+'_SINGS_wt'
fits_phot2 = path_cal+src+'_'+phot+'_SINGS'
tmp_phot2 = path_tmp+'calib/'+src+'_'+phot+'_SINGS'

fits_spec = path_cal+src+'_'+phot+'_IRC'
tmp_spec = path_tmp+'calib/'+src+'_'+phot+'_IRC'

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
                         dist='norm', Nmc=Nmc, filOUT=tmp_phot1+'_'+fname)

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
    
    ## (Sect. 3.2) The weight maps contain the information on the number of frames 
    ## that were used to create the science mosaics at each pixel (value= # frames x10); 
    ## the pixel size of the weight maps is the same as the science mosaics.
    ## -- SINGS v5 release doc
    ## https://irsa.ipac.caltech.edu/data/SPITZER/SINGS/doc/sings_fifth_delivery_v2.pdf
    
    ## Only need to run ONCE for each photometry (comment after using)
    # gen_unc2 = input("Create SINGS IRAC1 uncertainty map (y/n)? ")
    # if gen_unc2=='y':
    #     ds_wgt = read_fits(wgt_phot2)
    #     bg = read_fits(raw_phot2).data[694:734,56:96]
    #     bg_wgt = ds_wgt.data[694:734,56:96]
    #     wfac = .1
    #     iuncert(raw_phot2, filOUT=raw_phot2+'_unc',
    #             filWGT=wgt_phot2, wfac=wfac,
    #             BG_image=bg, BG_weight=bg_wgt)

    for f in fits_irc:
        ## Reproject phot
        fname = os.path.basename(f)
        refheader = fixwcs(f+fitsext).header
        mtg = imontage('exact', tmpdir=path_tmp)
        mtg.reproject_mc(raw_phot2, refheader=refheader,
                         dist='norm', Nmc=Nmc, filOUT=tmp_phot2+'_'+fname)

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

## Fit correlation
##-----------------
for f in fits_irc:
    fname = os.path.basename(f)

    ## Data
    data_phot1 = read_fits(fits_phot1+'_'+fname).data
    unc_phot1 = read_fits(fits_phot1+'_'+fname+'_unc').data
    data_phot2 = read_fits(fits_phot2+'_'+fname).data
    unc_phot2 = read_fits(fits_phot2+'_'+fname+'_unc').data
    data_spec = read_fits(fits_spec+'_'+fname).data
    unc_spec = read_fits(fits_spec+'_'+fname+'_unc').data
            
    if int_calib=='y':
        
        ## Plot
        pix_phot1 = data_phot1[:,:].reshape((-1,))
        pix_phot1_unc = unc_phot1[:,:].reshape((-1,))
        pix_phot2 = data_phot2[:,:].reshape((-1,))
        pix_phot2_unc = unc_phot2[:,:].reshape((-1,))
        pix_spec = data_spec[:,:].reshape((-1,))
        pix_spec_unc = unc_spec[:,:].reshape((-1,))
        
        # xgrid = np.arange(0,1e3,1)
        xgrid = np.logspace(-6,4,10000)
        
        ## NaN mask
        mask = ~np.logical_or( np.isnan(pix_spec), np.isnan(pix_phot2) )
        
        ## S/N ratio
        # print('S/N (IRC) = \n', pix_spec[mask]/pix_spec_unc[mask])
        # print('S/N (SINGS) = \n', pix_phot2[mask]/pix_phot2_unc[mask])
        # exit()

        if not '5124077' in fname:
            calib_gain = data_phot2/data_spec
            print(fname+' inter-calibration gain: {}'.format(calib_gain))
            
            ## Spectral correction
            ##---------------------
            if do_correct=='y':
                for j in trange(Nmc+1, leave=False,
                                desc=fname+': Spectral correction [MC]'):
                    if j==0:
                        imp = improve(f)
                    else:
                        imp = improve(f+'_'+str(j))
                    im = imp.im * calib_gain
                    ## Cover the original slits; Rebuild if need uncalibrated slits
                    write_fits(f+'_'+str(j), imp.hdr, im, imp.wvl, imp.wmod)
            else:
                sc = intercalib(f)
                sc.correct_spec(filOUT=f+'_0')
        else:
            sc = intercalib(f)
            sc.correct_spec(filOUT=f+'_0')

##----------------------------------------------------------

##                       Coadd slits

##----------------------------------------------------------
irc_coadd = input("Coadd IRC slits (y/n)? ")
if irc_coadd=='y':

    ## Coadding
    ##----------

    ## IRC grid
    refheader = fixwcs(fits_irc[0]+fitsext).header
    swp = iswarp(fits_irc, refheader=refheader, 
                 tmpdir=path_build, verbose=verbose)
    coadd_fp = swp.refheader

    if coadd_tool=='swarp':
        
        ## <iswarp> coadding
        ##===================
        swp = iswarp(refheader=coadd_fp,
                     # center='9:55:52,69:40:45', pixscale=6,
                     tmpdir=path_build, verbose=verbose)
    
        for j in trange(Nmc+1, leave=False, 
                        desc='<iswarp> IRC Coadding [MC]'):

            fits_irc_mc = []
            for f in fits_irc:
                fits_irc_mc.append(f+'_'+str(j))
            comb = swp.combine(fits_irc_mc,# combtype='wgt_avg',
                               keepedge=True, cropedge=True,
                               tmpdir=path_build+'MC_'+str(j)+'/',
                               filOUT=out_irc+'_'+str(j))
            ## Rescale pixel size
            if j==Nmc:
                rebinfo = True
            else:
                rebinfo = verbose
            irebin(out_irc+'_'+str(j), filOUT=out_irc+'_'+str(j),
                   pixscale=pixscale,
                   total=False, extrapol=True, verbose=rebinfo)


    elif coadd_tool=='reproject':

        ## <imontage> coadding
        ##=====================
        mtg = imontage('exact', tmpdir=path_tmp, verbose=verbose)

        for j in trange(Nmc+1, leave=False,
                        desc='<imontage> IRC coadding [MC]'):
            fits_irc_mc = []
            for f in fits_irc:
                fits_irc_mc.append(f+'_'+str(j))
            mtg.coadd(fits_irc_mc, refheader=coadd_fp,
                      filOUT=out_irc+'_'+str(j))
            ## Rescale pixel size
            if j==Nmc:
                rebinfo = True
            else:
                rebinfo = verbose
            irebin(out_irc+'_'+str(j), filOUT=out_irc+'_'+str(j),
                   pixscale=pixscale,
                   total=False, extrapol=True, verbose=rebinfo)
        
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
        
##----------------------------------------------------------

##                 Atlas inter-Calibration

##----------------------------------------------------------
pre_calib1_all = input("Run atlas inter-calibration prepipeline (DustPedia) (y/n)? ")
pre_calib2_all = input("Run atlas inter-calibration prepipeline (SINGS) (y/n)? ")
int_calib_all = input("Run atlas inter-calibration (y/n)? ")

## Prepare photometry
##--------------------

## DustPedia (phot1)
##===================
if pre_calib1_all=='y':
    
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
if pre_calib2_all=='y':
    
    ## Create uncertainty via weight map (suppose uniform contribution)
    
    ## (Sect. 3.2) The weight maps contain the information on the number of frames 
    ## that were used to create the science mosaics at each pixel (value= # frames x10); 
    ## the pixel size of the weight maps is the same as the science mosaics.
    ## -- SINGS v5 release doc
    ## https://irsa.ipac.caltech.edu/data/SPITZER/SINGS/doc/sings_fifth_delivery_v2.pdf
    
    ## Only need to run ONCE for each photometry (comment after using)
    # gen_unc2 = input("Create SINGS IRAC1 uncertainty map (y/n)? ")
    # if gen_unc2=='y':
    #     ds_wgt = read_fits(wgt_phot2)
    #     bg = read_fits(raw_phot2).data[694:734,56:96]
    #     bg_wgt = ds_wgt.data[694:734,56:96]
    #     wfac = .1
    #     iuncert(raw_phot2, filOUT=raw_phot2+'_unc',
    #             filWGT=wgt_phot2, wfac=wfac,
    #             BG_image=bg, BG_weight=bg_wgt)

    ## Reproject phot
    refheader = fixwcs(out_irc+'_0'+fitsext).header
    mtg = imontage('exact', tmpdir=path_tmp)
    mtg.reproject_mc(raw_phot2, refheader=refheader,
                     dist='norm', Nmc=Nmc, filOUT=tmp_phot2)

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

## Fit correlation
##-----------------
## Data
data_phot1 = read_fits(fits_phot1).data
unc_phot1 = read_fits(fits_phot1+'_unc').data
data_phot2 = read_fits(fits_phot2).data
unc_phot2 = read_fits(fits_phot2+'_unc').data
data_spec = read_fits(fits_spec).data
unc_spec = read_fits(fits_spec+'_unc').data

if int_calib_all=='y':
    
    ## Plot
    pix_phot1 = data_phot1[:,:].reshape((-1,))
    pix_phot1_unc = unc_phot1[:,:].reshape((-1,))
    pix_phot2 = data_phot2[:,:].reshape((-1,))
    pix_phot2_unc = unc_phot2[:,:].reshape((-1,))
    pix_spec = data_spec[:,:].reshape((-1,))
    pix_spec_unc = unc_spec[:,:].reshape((-1,))
    
    # xgrid = np.arange(0,1e3,1)
    xgrid = np.logspace(-6,4,10000)
    
    ## NaNs mask
    mask = ~np.logical_or( np.isnan(pix_spec), np.isnan(pix_phot2) )

    ## SINGS - DP
    p0 = pplot(pix_phot1, pix_phot2,
               yerr=pix_phot2_unc, xerr=pix_phot1_unc,
               fmt='.', c='k', ec='r', elw=1,
               xlog=1, ylog=1, nonposx='clip', nonposy='clip',
               xlim=(1e-3,1e3), ylim=(1e-3,1e3),
               xlabel='DustPedia (MJy/sr)', ylabel='SINGS (MJy/sr)',
               figsize=(8,8), legend='upper left', title=src+'_'+phot,
               titlesize=20, labelsize=10, ticksize=10, legendsize=10)
    p0.save(path_cal+'SINGS-DP_'+phot+'.png')
            
    ## IRC - SINGS
    p = pplot(pix_phot2, pix_spec,
              yerr=pix_spec_unc, xerr=pix_phot2_unc,
              fmt='.', c='k', ec='r', 
              xlog=1, ylog=1, nonposx='clip', nonposy='clip',
              xlim=(1e-6,1e3), ylim=(1e-6,1e3),
              xlabel='SINGS (MJy/sr)', ylabel='IRC (MJy/sr)', 
              figsize=(8,8), legend='upper left', title=src+'_'+phot,
              titlesize=20, labelsize=10, ticksize=10, legendsize=10)

    ## Linear fit
    popt, pcov = curve_fit(f_lin0, pix_spec[mask], pix_phot2[mask],
                                   sigma=pix_phot2_unc[mask])
    calib_gain = popt[0]
    # calib_off = popt[1]
    calib_off = 0.
    # calib_gain = 1.
    # calib_off = popt[0]
    print('Inter-Calibration ('+phot+') gain = {:.4}'.format(calib_gain))
    print('Inter-Calibration ('+phot+') offset = {:.4}'.format(calib_off))
    label = phot+': y={0:.4}x+{1:.4}'.format(calib_gain, calib_off)
    p.add_plot(f_lin0(xgrid, *popt), xgrid,
               c='y', ls='-', label=label)
    
    # popt, pcov = curve_fit(f_lin0, pix_phot2[mask], pix_spec[mask],)
    #                        # sigma=pix_spec_unc[mask])
    # calib_gain = 1. / popt[0]
    # # calib_off = -popt[1] / popt[0]
    # calib_off = 0.
    # print('Inter-Calibration ('+phot+') gain = {:.4}'.format(calib_gain))
    # print('Inter-Calibration ('+phot+') offset = {:.4}'.format(calib_off))
    # label = phot+': y={0:.4}x+{1:.4}'.format(calib_gain, calib_off)
    # p.add_plot(xgrid, f_lin0(xgrid, *popt),
    #            c='y', ls='-', label=label)
    
    p.save(path_cal+'intercalib_'+phot+'.png')
            
    ## Spectral correction
    ##---------------------    
    sc = intercalib(out_irc+'_0')
    sc.correct_spec(filOUT=out_irc) # Duplicate spectral cube
    sc.read_filter(phot) # Convert calib_off to spec_off
    spec_off = calib_off * sc.specoff_ov_bboff[0]
    spec_gain = calib_gain
        
    do_correct_all = input(" - Do IRC spectral correction [{}] (y/n)? ".format(phot))
    if do_correct_all=='y':
        do_indiv = input("   - Do individual correction for each pixel (y/n)? ")

        if do_indiv=='y':
            ## Individual
            ##============
            spec_gain = (data_phot2-spec_off)/data_spec
            mask2D = ~np.logical_or( np.isnan(data_spec), np.isnan(data_phot2) )
            imp = improve(out_irc)
            im = imp.im
            for k, lam in enumerate(imp.wvl):
                im[k][mask2D] = imp.im[k][mask2D] * spec_gain[mask2D] + spec_off
                
            write_fits(out_irc, imp.hdr, im, imp.wvl, imp.wmod)
            
        else:
            ## Global
            ##========
            sc = intercalib(out_irc)
            sc.correct_spec(filOUT=out_irc,
                            gain=spec_gain, offset=spec_off)
        
##----------------------------------------------------------

##                      Plot spectra

##----------------------------------------------------------
do_plot = input("Plot IRC spectra (y/n)? ")
if do_plot=='y':
    ds = read_fits(out_irc, out_irc+'_unc')
    # ds.data = respect().smooth(out_irc+'_0', out_irc+'_0', lim_unc=1.e2, cmin=2)
    ds0 = read_fits(out_irc+'_0')
    pp = intercalib(out_irc)
    pp.read_filter(phot)
    wcen = pp.wcen
    Nw, Ny, Nx = ds.data.shape
    for x in range(Nx):
        for y in range(Ny):
            if ~np.isnan(ds.data[:,y,x]).any():
                p = pplot(ds.wave, ds.data[:,y,x], yerr=ds.unc[:,y,x],
                          # xlog=1, ylog=1, 
                          c='k', lw=.7, ec='r', 
                          figsize=(8,8), legend='upper left',
                          title='IRC_('+str(x+1)+','+str(y+1)+')',
                          titlesize=20, labelsize=10, ticksize=10, legendsize=10)

                p.add_plot(ds0.wave, ds0.data[:,y,x],
                           c='y', lw=.7, ls='--', zorder=-1)
                
                p.add_plot(wcen[0], data_phot2[y,x], yerr=unc_phot2[y,x],
                           c='m', marker='o', ms=10, zorder=100, label=phot)
                p.add_plot(wcen[0], data_spec[y,x], yerr=unc_spec[y,x],
                           c='g', marker='^', ms=10, zorder=101, label='IRC-'+phot)
                
                p.save(path_fig+'IRC_('+str(x+1)+','+str(y+1)+')')

