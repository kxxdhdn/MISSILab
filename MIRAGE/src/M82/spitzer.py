#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""

Build stitched Spitzer/IRS spectral cube
Inter-calibrate Spitzer/IRS spectra with Spitzer/IRAC4 (SL) & MIPS1 (LL, SH+LH)
Combine AKARI IRC slits

"""

# import logging, sys
# logging.basicConfig(stream=sys.stderr, level=logging.ERROR)
# print(logging.getLogger())
# logging.disable(sys.maxsize) # disable IDL print
import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning) 

from tqdm import tqdm, trange

import os, pathlib
import math
import numpy as np
from scipy.optimize import curve_fit
from scipy.integrate import trapz

## rapyuta
from rapyuta.inout import fclean, fitsext, read_fits, write_fits, read_hdf5
from rapyuta.imaging import ( iconvolve, iswarp, iuncert, icrop,
                              igroupixel, concatenate, imontage,
                              improve, Jy_per_pix_to_MJy_per_sr )
from rapyuta.astrom import fixwcs
from rapyuta.arrays import closest, pix2sup, sup2pix
from rapyuta.calib import intercalib
from rapyuta.maths import f_lin, f_lin0, f_lin1
from rapyuta.plots import pplot


##----------------------------------------------------------

##                      Preliminaries

##----------------------------------------------------------

## Local
from buildinfo import ( src, Nmc, verbose, coadd_tool, colors, markers,
                        chnl, path_irs, sub_SL, sub_LL, sub_SH, sub_LH,
                        path_idl, path_ker, path_conv, path_phot, path_cal,
                        fits_ker, fwhm, csv_ker, path_tmp, path_fig, path_out,
                        filog, fits_irc, out_irc, out_irs, path_build, parobs, # build
)
from wcsinfo import coadd_footprint


## Banner
print('\n============================================================\n')

print('        MIRAGE - Spitzer/IRS cube builder - '+src)

print('\n============================================================\n')

# Nch = len(chnl)
chnl_SL = []
for c in chnl:
    if c[:2]=='SL':
        chnl_SL.append(c)
Nch_SL = len(chnl_SL)
chnl_LL = []
for c in chnl:
    if c[:2]=='LL':
        chnl_LL.append(c)
Nch_LL = len(chnl_LL)

Nsub_SL = len(sub_SL)
Nsub_LL = len(sub_LL)
Nsub_SH = len(sub_SH)
Nsub_LH = len(sub_LH)

if ('SH' in chnl) or ('LH' in chnl):
    phot = ['MIPS1']
else:
    phot = ['IRAC4', 'MIPS1']
Nphot = len(phot)

Nobs = len(out_irc)

# Nmc = 2

resume = False

# exit()


##----------------------------------------------------------

##              Sub-map processing & Stitching

##----------------------------------------------------------
proc_SL = input("Process SL (y/n)? ")
proc_LL = input("Process LL (y/n)? ")
do_stitch = input("Stitch SL/LL (y/n) ")
if do_stitch=='y':
    match_SL = input("Match SL2-SL1 (y/n)? ")
    match_LL = input("Match LL2-LL1 (y/n)? ")
# proc_SH = input("Process SH (y/n)? ")
# proc_LH = input("Process LH (y/n)? ")
# match_hires = input(" - Match SH-LH (y/n)? ")
proc_SH = None
proc_LH = None
match_hires = None

## Prepare headers of sub-maps
##-----------------------------
if proc_SL=='y':
    hdr_SL = []
    for i in trange(Nsub_SL,# leave=False,
                    desc='Loading SL refheader'):
        os.makedirs(path_tmp+'SL/', exist_ok=True)
        ref0 = []
        for c in chnl_SL:
            ds = read_fits(path_irs+src+'_'+sub_SL[i]+'_'+c)
            write_fits(path_tmp+'SL/tmpref_'+sub_SL[i]+'_'+c, ds.header, ds.data[:2], ds.wave[:2])
            ref0.append(path_tmp+'SL/tmpref_'+sub_SL[i]+'_'+c)
        swp0 = iswarp(ref0, tmpdir=path_tmp+'SL/')
        ref = []
        for c in chnl_SL:
            swp0.combine(path_tmp+'SL/tmpref_'+sub_SL[i]+'_'+c,
                         keepedge=True,# cropedge=True,
                         filOUT=path_tmp+'SL/tmpref_'+c)
            ref.append(path_tmp+'SL/tmpref_'+c)
        hdr = concatenate(ref, cropedge=True,
                          wsort=False, keepfrag=False).header
        hdr_SL.append(fixwcs(header=hdr).header)
    fclean(path_tmp+'SL/tmpref_*')
if proc_LL=='y':
    hdr_LL = []
    for i in trange(Nsub_LL,# leave=False,
                    desc='Loading LL refheader'):
        os.makedirs(path_tmp+'LL/', exist_ok=True)
        ref0 = []
        for c in chnl_LL:
            ds = read_fits(path_irs+src+'_'+sub_LL[i]+'_'+c)
            write_fits(path_tmp+'LL/tmpref_'+sub_LL[i]+'_'+c, ds.header, ds.data[:2], ds.wave[:2])
            ref0.append(path_tmp+'LL/tmpref_'+sub_LL[i]+'_'+c)
        swp0 = iswarp(ref0, tmpdir=path_tmp+'LL/')
        ref = []
        for c in chnl_LL:
            swp0.combine(path_tmp+'LL/tmpref_'+sub_LL[i]+'_'+c,
                         keepedge=True,# cropedge=True,
                         filOUT=path_tmp+'LL/tmpref_'+c)
            ref.append(path_tmp+'LL/tmpref_'+c)
        hdr = concatenate(ref, cropedge=True,
                          wsort=False, keepfrag=False).header
        hdr_LL.append(fixwcs(header=hdr).header)
    fclean(path_tmp+'LL/tmpref_*')

## Coadd SL/LL modules (using the bonus channel)
##-----------------------------------------------
if do_stitch=='y':
    for iph in range(Nphot):
        if phot[iph]=='IRAC4':
            sub_XX = sub_SL
            path_tmp_XX = path_tmp+'SL/'
            mod = 'SL'
        elif phot[iph]=='MIPS1':
            sub_XX = sub_LL
            path_tmp_XX = path_tmp+'LL/'
            mod = 'LL'
        os.makedirs(path_tmp_XX, exist_ok=True)
    
        if resume:
            iresume = int(input("Resume "+mod+" sub-map processing from the iteration: "))
        else:
            iresume = 0
        for j in trange(Nmc+1-iresume,# leave=False,
                        desc='IRS '+mod+' sub-map processing [MC]'):
            j += iresume
            
            for i in trange(len(sub_XX), leave=False,
                            desc=' - Sub-maps'):
                if np.logical_or(proc_SL=='y' and phot[iph]=='IRAC4',
                                 proc_LL=='y' and phot[iph]=='MIPS1'):
                    if phot[iph]=='IRAC4':
                        refheader = hdr_SL[i]
                        chnl_XX = chnl_SL
                    elif phot[iph]=='MIPS1':
                        refheader = hdr_LL[i]
                        chnl_XX = chnl_LL
                        
                    for ich in trange(len(chnl_XX), leave=False,
                                      desc='   - '+src+'_'+sub_XX[i]+'_'+mod):
                        f0nam = src+'_'+sub_XX[i]+'_'+chnl_XX[ich]
                        f0nam_in = path_irs+f0nam
                        f0nam_out =  path_tmp_XX+f0nam+'_'+str(j)
        
                        ## Reprojection 
                        ##--------------
                        if coadd_tool=='swarp':
                            swp = iswarp(refheader=refheader, tmpdir=path_tmp_XX)
                            if j==0:
                                swp.combine(f0nam_in, keepedge=True,# cropedge=True,
                                            filOUT=f0nam_out)
                            else:
                                swp.combine(f0nam_in, keepedge=True,# cropedge=True,
                                            dist='norm', sig_pt=.2,
                                            filOUT=f0nam_out)
                        elif coadd_tool=='reproject':
                            mtg = imontage('exact', tmpdir=path_tmp_XX)
                            if j==0:
                                mtg.coadd(f0nam_in, refheader=refheader,
                                          filOUT=f0nam_out)
                            else:
                                mtg.coadd(f0nam_in, refheader=refheader,
                                          dist='norm', sig_pt=.2,
                                          filOUT=f0nam_out)
        
                        ## PSF Convolution with both sources of errors (MC)
                        ##--------------------------------------------------
                        conv = iconvolve(f0nam_out, kfile=fits_ker, psf=fwhm,
                                         klist=csv_ker, convdir=path_conv,
                                         filOUT=f0nam_out)
                        conv.do_conv(idldir=path_idl)
        
                ##-----------------------------
                ## SL2-SL1 / LL2-LL1 stitching
                ##-----------------------------
                if phot[iph]=='IRAC4':
                    f1nam = path_tmp+'SL/'+src+'_'+sub_SL[i]+'_SL'
                elif phot[iph]=='MIPS1':
                    f1nam = path_tmp+'LL/'+src+'_'+sub_LL[i]+'_LL'
                
                if np.logical_or(match_SL=='y' and phot[iph]=='IRAC4',
                                 match_LL=='y' and phot[iph]=='MIPS1'):
                    ## Read spectra
                    data3 = read_fits(f1nam+'3_'+str(j)).data
                    data2 = read_fits(f1nam+'2_'+str(j)).data
                    data1 = read_fits(f1nam+'1_'+str(j)).data
                    wvl3 = read_fits(f1nam+'3_'+str(j)).wave
                    wvl2 = read_fits(f1nam+'2_'+str(j)).wave
                    wvl1 = read_fits(f1nam+'1_'+str(j)).wave
                    
                    ## Match SL2 to SL3 (idem. LL)
                    ## --iwmin_SL3--iwmin_SL1--iwmax_SL2--iwmax_SL3--
                    ##       |          |          |          |
                    ##       |------left_SL3-------|
                    ##       |------right_SL2------|
                    iwmax2 = closest(wvl3, wvl2[-1], side='left') + 1 # SL3 index
                    left3 = trapz(data3[:iwmax2], wvl3[:iwmax2],
                                     dx=wvl3[1]-wvl3[0], axis=0)
                    iwmin3 = closest(wvl2, wvl3[0], side='right') # SL2 index
                    right2 = trapz(data2[iwmin3:], wvl2[iwmin3:],
                                      dx=wvl2[1]-wvl2[0], axis=0)
                    
                    ## Match SL1 to SL3 (idem. LL)
                    ## --iwmin_SL3--iwmin_SL1--iwmax_SL2--iwmax_SL3--
                    ##       |          |          |          |
                    ##                  |------right_SL3------|
                    ##                  |------left_SL1-------|
                    iwmin1 = closest(wvl3, wvl1[0], side='right') # SL3 index
                    right3 = trapz(data3[iwmin1:], wvl3[iwmin1:],
                                      dx=wvl3[1]-wvl3[0], axis=0)
                    iwmax3 = closest(wvl1, wvl3[-1], side='left') + 1 # SL1 index
                    left1 = trapz(data1[:iwmax3], wvl1[:iwmax3],
                                     dx=wvl1[1]-wvl1[0], axis=0)
        
                    ## Calculate scaling factors
                    gain2 = left3 / right2
                    gain1 = right3 / left1
                    ## Display scaling factor map
                    mask2D = ~np.isnan(gain2)
                    # print('SL2 to SL3 scaling factor: ', gain2[mask2D])
                    mask2D = ~np.isnan(gain1)
                    # print('SL1 to SL3 scaling factor: ', gain1[mask2D])
                else:
                    gain2 = 1
                    gain1 = 1
    
                fits_ord2 = pathlib.Path(f1nam+'2_'+str(j)+fitsext)
                fits_ord1 = pathlib.Path(f1nam+'1_'+str(j)+fitsext)
                if fits_ord2.exists() and fits_ord1.exists():
                    ic = intercalib(f1nam+'2_'+str(j))
                    ic.correct_spec(gain=gain2, filOUT=f1nam+'2_match')
                    ic = intercalib(f1nam+'1_'+str(j))
                    ic.correct_spec(gain=gain1, filOUT=f1nam+'1_match')
                    ## Stitch
                    concatenate( (f1nam+'2_match',f1nam+'1_match'),
                                 f1nam+'_'+str(j), wsort=False,
                                 keepfrag=False, cropedge=False )
                else:
                    if fits_ord2.exists():
                        ic = intercalib(f1nam+'2_'+str(j))
                        ic.correct_spec(gain=gain2, filOUT=f1nam+'_'+str(j))
                    if fits_ord1.exists():
                        ic = intercalib(f1nam+'1_'+str(j))
                        ic.correct_spec(gain=gain1, filOUT=f1nam+'_'+str(j))

## Coadd hires modules (using one-point/integral scaling factor)
##---------------------------------------------------------------
# if resume:
#     iresume = int(input("Resume sub-map processing from the iteration: "))
# else:
#     iresume = 0
for j in trange(Nmc+1,#-iresume,# leave=False,
                desc='IRS hires sub-map processing [MC]'):
    # j += iresume
    
    ##====
    ## SH
    ##====
    if proc_SH=='y':
        for i in trange(Nsub_SH, leave=False,
                        desc=' - SH sub-maps'):
            f0nam = src+'_'+sub_SH[i]+'_SH'

            ## PSF Convolution with both sources of errors (MC)
            ##--------------------------------------------------
            if j==0:
                conv = iconvolve(path_irs+f0nam, kfile=fits_ker, psf=fwhm,
                                 klist=csv_ker, convdir=path_conv,
                                 filOUT=path_tmp+'SH/'+f0nam+'_'+str(j))
            else:
                conv = iconvolve(path_irs+f0nam, kfile=fits_ker, psf=fwhm,
                                 klist=csv_ker, convdir=path_conv,
                                 dist='norm', sig_pt=.2,
                                 filOUT=path_tmp+'SH/'+f0nam+'_'+str(j))
            conv.do_conv(idldir=path_idl)
    
    ##====
    ## LH
    ##====
    if proc_LH=='y':
        for i in trange(Nsub_LH, leave=False,
                        desc=' - LH sub-maps'):
            f0nam = src+'_'+sub_LH[i]+'_LH'

            ## PSF Convolution with both sources of errors (MC)
            ##--------------------------------------------------
            if j==0:
                conv = iconvolve(path_irs+f0nam, kfile=fits_ker, psf=fwhm,
                                 klist=csv_ker, convdir=path_conv,
                                 filOUT=path_tmp+'LH/'+f0nam+'_'+str(j))
            else:
                conv = iconvolve(path_irs+f0nam, kfile=fits_ker, psf=fwhm,
                                 klist=csv_ker, convdir=path_conv,
                                 dist='norm', sig_pt=.2,
                                 filOUT=path_tmp+'LH/'+f0nam+'_'+str(j))
            conv.do_conv(idldir=path_idl)
            
    ##--------------
    ## Stitch SH-LH
    ##--------------
    if match_hires=='y':
        
        ## Coadd all SH sub-maps
        ##-----------------------
        fsub_SH = []
        for i in range(Nsub_SH):
            f1nam = path_tmp+'SH/'+src+'_'+sub_SH[i]+'_SH'
            fsub_SH.append(f1nam+'_'+str(j))
        refheader = iswarp(fsub_SH, tmpdir=path_tmp).refheader
        
        if coadd_tool=='swarp':
            swp = iswarp(refheader=refheader,
                         tmpdir=path_tmp, verbose=verbose)
            swp.combine(fsub_SH,# combtype='wgt_avg',
                        keepedge=True,# cropedge=True,
                        filOUT=path_tmp+src+'_SH_'+str(j))
        elif coadd_tool=='reproject':
            mtg = imontage('exact', tmpdir=path_tmp)
            mtg.coadd(fsub_SH, refheader=refheader,
                      filOUT=path_tmp+src+'_SH_'+str(j))

        ## Coadd all LH sub-maps
        ##-----------------------
        fsub_LH = []
        for i in range(Nsub_LH):
            f1nam = path_tmp+'LH/'+src+'_'+sub_LH[i]+'_LH'
            fsub_LH.append(f1nam+'_'+str(j))
        refheader = iswarp(fsub_LH, tmpdir=path_tmp).refheader
        
        if coadd_tool=='swarp':
            swp = iswarp(refheader=refheader,
                         tmpdir=path_tmp, verbose=verbose)
            swp.combine(fsub_LH,# combtype='wgt_avg',
                        keepedge=True,# cropedge=True,
                        filOUT=path_tmp+src+'_LH_'+str(j))
        elif coadd_tool=='reproject':
            mtg = imontage('exact', tmpdir=path_tmp)
            mtg.coadd(fsub_LH, refheader=refheader,
                      filOUT=path_tmp+src+'_LH_'+str(j))
        
        ## Coadd SH & LH maps
        ##---------------------
        refheader = iswarp((path_tmp+src+'_SH',path_tmp+src+'_LH'),
                           tmpdir=path_tmp+'LH/').refheader

        if coadd_tool=='swarp':
            swp = iswarp(refheader=refheader,
                         tmpdir=path_tmp, verbose=verbose)
            swp.combine(path_tmp+src+'_SH_'+str(j),# combtype='wgt_avg',
                        keepedge=True,# cropedge=True,
                        filOUT=path_tmp+src+'_SH_'+str(j))
            swp.combine(path_tmp+src+'_LH_'+str(j),# combtype='wgt_avg',
                        keepedge=True,# cropedge=True,
                        filOUT=path_tmp+src+'_LH_'+str(j))
        elif coadd_tool=='reproject':
            mtg = imontage('exact', tmpdir=path_tmp)
            mtg.coadd(path_tmp+src+'_SH_'+str(j), refheader=refheader,
                      filOUT=path_tmp+src+'_SH_'+str(j))
            mtg.coadd(path_tmp+src+'_LH_'+str(j), refheader=refheader,
                      filOUT=path_tmp+src+'_LH_'+str(j))
        
        ## SH
        data_SH = read_fits(path_tmp+src+'_SH_'+str(j)).data
        wvl_SH = read_fits(path_tmp+src+'_SH_'+str(j)).wave
        ## LH
        data_LH = read_fits(path_tmp+src+'_LH_'+str(j)).data
        wvl_LH = read_fits(path_tmp+src+'_LH_'+str(j)).wave
        
        iwmax_SH = closest(wvl_LH, wvl_SH[-1], side='left') + 1 # LH index
        iwmin_LH = closest(wvl_SH, wvl_LH[0], side='right') # SH index
        right_SH = trapz(data_SH[iwmin_LH:], wvl_SH[iwmin_LH:],
                         dx=wvl_SH[1]-wvl_SH[0], axis=0)
        left_LH = trapz(data_LH[:iwmax_SH], wvl_LH[:iwmax_SH],
                        dx=wvl_LH[1]-wvl_LH[0], axis=0)
    
        ## Integral scaling factor
        gain_hires = right_SH / left_LH
        
        ## One-point scaling factor
        # i_SH = closest(wvl_SH, 19.10)
        # i_LH = closest(wvl_LH, 19.30)
        # gain_hires = data_SH[i_SH] / data_LH[i_LH]
        
        ## Display scaling factor map
        mask2D = ~np.isnan(gain_hires)
        # print('LH to SH scaling factor: ', gain_hires[mask2D])

        ic = intercalib(path_tmp+src+'_LH_'+str(j))
        ic.correct_spec(gain=gain_hires,
                        filOUT=path_tmp+src+'_LH_match')

        wran_hi = [ (10.00, 19.19), # SH
                    (19.20, 37.10), ] # LH
        concatenate((path_tmp+src+'_SH_'+str(j),path_tmp+src+'_LH_match'),
                    path_tmp+src+'_hires'+str(j), wsort=False, wrange=wran_hi,
                    keepfrag=False, cropedge=False)

if (proc_SH=='y' or proc_SH=='y' or match_hires=='y'):
    print('Plot spectra...')
    exit()

##----------------------------------------------------------

##          Inter-Calibration (sub-maps) & Coadding

##----------------------------------------------------------
## DustPedia has no MIPS1 available for M82
## SINGS MIPS1 map of M82 has a saturation hole in the central disk
photcal = input("Photometry choice: SINGS/DustPedia (s/d)? ")
pre_calib = input("Prepare photometry (y/n)? ")
synt_phot = input("Prepare IRS synthetic photometry (y/n)? ")
write_SL = input("Write IRS-SL spectral sub-maps (y/n)? ")
coadd_SL = input("Coadd IRS-SL spectral sub-maps (y/n)? ")
build_SL = input("Build IRS-SL slits (y/n)? ")
write_LL = input("Write IRS-LL spectral sub-maps (y/n)? ")
coadd_LL = input("Coadd IRS-LL spectral sub-maps (y/n)? ")
build_LL = input("Build IRS-LL slits (y/n)? ")

os.makedirs(path_tmp+'calib/', exist_ok=True)

## Sub-map homogenize photometry and synthetic photometry
##--------------------------------------------------------
for iph in range(Nphot):

    ## Prepare photometry
    if phot[iph]=='IRAC4':
        phot_ker = path_ker+'Kernel_HiRes_IRAC_8.0_to_Gauss_06.0'
    elif phot[iph]=='MIPS1':
        phot_ker = path_ker+'Kernel_HiRes_MIPS_24_to_Gauss_06.0'

    if photcal=='d' and phot[iph]=='IRAC4':
        cal_phot = path_cal+src+'_'+phot[iph]+'_DP'
        tmp_phot = path_tmp+'calib/'+src+'_'+phot[iph]+'_DP'
    else:
        cal_phot = path_cal+src+'_'+phot[iph]+'_SINGS'
        tmp_phot = path_tmp+'calib/'+src+'_'+phot[iph]+'_SINGS'
    cal_spec = path_cal+src+'_'+phot[iph]+'_IRS'
    tmp_spec = path_tmp+'calib/'+src+'_'+phot[iph]+'_IRS'
    
    ##------------------------
    ## SINGS/DustPedia (phot)
    ##------------------------
    if pre_calib=='y':
        if resume:
            iresume = int(input("Resume inter-calib prepipeline ("+phot[iph]+
                                ") from the iteration: "))
        else:
            iresume = 0

        if photcal=='d' and phot[iph]=='IRAC4':
            raw_phot = path_phot+src+'_'+phot[iph]+'_DP'

            ## Convert phot unit (Jy/pix -> MJy/sr)
            ## This step should be before reprojection (pixscale may change)
            Jy_per_pix_to_MJy_per_sr(raw_phot, filOUT=tmp_phot)

        else:
            raw_phot = path_phot+src+'_'+phot[iph]+'_SINGS'

            ## Create uncertainty via weight map (suppose uniform contribution)
            
            ## (Sect. 3.2) The weight maps contain the information on the number of frames
            ## that were used to create the science mosaics at each pixel (value= # frames x10);
            ## the pixel size of the weight maps is the same as the science mosaics.
            ## -- SINGS v5 release doc
            ## https://irsa.ipac.caltech.edu/data/SPITZER/SINGS/doc/sings_fifth_delivery_v2.pdf
            
            ## Only need to run ONCE for each photometry (comment after using)
            # gen_unc = input("Create SINGS IRAC1 uncertainty map (y/n)? ")
            # wgt_phot = path_phot+src+'_'+phot[iph]+'_SINGS_wt'
            # if gen_unc=='y':
            #     ds_wgt = read_fits(wgt_phot)
            #     bg = read_fits(raw_phot).data[694:734,56:96]
            #     bg_wgt = ds_wgt.data[694:734,56:96]
            #     wfac = .1
            #     iuncert(raw_phot, filOUT=raw_phot+'_unc',
            #             filWGT=wgt_phot, wfac=wfac,
            #             BG_image=bg, BG_weight=bg_wgt)
        
        if phot[iph]=='IRAC4' and Nch_SL>1:
            ##====
            ## SL
            ##====
            for i in trange(Nsub_SL,# leave=False,
                            desc=phot[iph]+' processing'):
                f1nam = path_tmp+'SL/'+src+'_'+sub_SL[i]+'_SL'
                p1nam0 = cal_phot+'_'+sub_SL[i]+'_SL'
                refheader = fixwcs(f1nam+'_0'+fitsext).header

                for j in trange(Nmc+1-iresume, leave=False,
                                desc=' - '+src+'_'+sub_SL[i]+'_SL [MC]'):
                    j += iresume
                    
                    p1namj = tmp_phot+'_'+sub_SL[i]+'_SL_'+str(j)
                    
                    ## Reproject phot
                    if coadd_tool=='swarp':
                        swp = iswarp(refheader=refheader, tmpdir=path_tmp)
                        if j==0:
                            swp.combine(raw_phot, keepedge=True,
                                        filOUT=p1nam0)
                        else:
                            swp.combine(raw_phot, dist='norm', keepedge=True,
                                        filOUT=p1namj)
                    elif coadd_tool=='reproject':
                        mtg = imontage('exact', tmpdir=path_tmp)
                        if j==0:
                            mtg.reproject(raw_phot, refheader=refheader,
                                          filOUT=p1nam0)
                        else:
                            mtg.reproject(raw_phot, refheader=refheader, dist='norm',
                                          filOUT=p1namj)
                            
                    ## Convolve phot
                    if j==0:
                        conv = iconvolve(p1nam0,
                                         kfile=phot_ker, klist=csv_ker,
                                         filOUT=p1nam0)
                    else:
                        conv = iconvolve(p1namj,
                                         kfile=phot_ker, klist=csv_ker,
                                         filOUT=p1namj)
                    conv.do_conv(idldir=path_idl)
                    
        elif phot[iph]=='MIPS1' and Nch_LL>1:
            ##====
            ## LL
            ##====
            for i in trange(Nsub_LL,# leave=False,
                            desc=phot[iph]+' (SINGS) processing'):
                f1nam = path_tmp+'LL/'+src+'_'+sub_LL[i]+'_LL'
                p1nam0 = cal_phot+'_'+sub_LL[i]+'_LL'
                refheader = fixwcs(f1nam+'_0'+fitsext).header

                for j in trange(Nmc+1-iresume, leave=False,
                                desc=' - '+src+'_'+sub_LL[i]+'_LL [MC]'):
                    j += iresume
                    
                    p1namj = tmp_phot+'_'+sub_LL[i]+'_LL_'+str(j)
                    
                    ## Reproject phot
                    if coadd_tool=='swarp':
                        swp = iswarp(refheader=refheader, tmpdir=path_tmp)
                        if j==0:
                            swp.combine(raw_phot, keepedge=True,
                                        filOUT=p1nam0)
                        else:
                            swp.combine(raw_phot, dist='norm', keepedge=True,
                                        filOUT=p1namj)
                    elif coadd_tool=='reproject':
                        mtg = imontage('exact', tmpdir=path_tmp)
                        if j==0:
                            mtg.reproject(raw_phot, refheader=refheader,
                                          filOUT=p1nam0)
                        else:
                            mtg.reproject(raw_phot, refheader=refheader, dist='norm',
                                          filOUT=p1namj)
                            
                    ## Convolve phot
                    if j==0:
                        conv = iconvolve(p1nam0,
                                         kfile=phot_ker, klist=csv_ker,
                                         filOUT=p1nam0)
                    else:
                        conv = iconvolve(p1namj,
                                         kfile=phot_ker, klist=csv_ker,
                                         filOUT=p1namj)
                    conv.do_conv(idldir=path_idl)

    for i in range(Nsub_SL):
        p1nam0 = cal_phot+'_'+sub_SL[i]+'_SL'
        ## Calculate unc (SL)
        fits_Nmc = pathlib.Path(tmp_phot+'_'+sub_SL[i]+'_SL_'+str(Nmc)+fitsext)
        if Nmc>1 and fits_Nmc.exists():
            mcimage = []
            for j in range(Nmc):
                ds = read_fits(tmp_phot+'_'+sub_SL[i]+'_SL_'+str(j+1))
                mcimage.append(ds.data)
            mcimage = np.array(mcimage)
            unc = np.nanstd(mcimage, axis=0)
            write_fits(p1nam0+'_unc', ds.header, unc)

    for i in range(Nsub_LL):
        p1nam0 = cal_phot+'_'+sub_LL[i]+'_LL'
        ## Calculate unc (LL)
        fits_Nmc = pathlib.Path(tmp_phot+'_'+sub_LL[i]+'_LL_'+str(Nmc)+fitsext)
        if Nmc>1 and fits_Nmc.exists():
            mcimage = []
            for j in range(Nmc):
                ds = read_fits(tmp_phot+'_'+sub_LL[i]+'_LL_'+str(j+1))
                mcimage.append(ds.data)
            mcimage = np.array(mcimage)
            unc = np.nanstd(mcimage, axis=0)
            write_fits(p1nam0+'_unc', ds.header, unc)

    ##-----------------------------
    ## Synthetic photometry (spec)
    ##-----------------------------
    if synt_phot=='y':
        if resume:
            iresume = int(input("Resume synthetic photometry "+phot[iph]+
                                " from the iteration: "))
        else:
            iresume = 0
        
        if phot[iph]=='IRAC4' and Nch_SL>1:
            ##====
            ## SL
            ##====
            for i in trange(Nsub_SL,# leave=False,
                            desc=phot[iph]+' (SL synt phot) sub-maps'):
                f1nam = path_tmp+'SL/'+src+'_'+sub_SL[i]+'_SL'
                s1nam0 = cal_spec+'_'+sub_SL[i]+'_SL'

                for j in trange(Nmc+1-iresume, leave=False,
                                desc=' - '+sub_SL[i]+'_SL [MC]'):
                    j += iresume

                    s1namj = tmp_spec+'_'+sub_SL[i]+'_SL_'+str(j)
                    
                    ic = intercalib(f1nam+'_'+str(j))
                    sp = ic.synthetic_photometry(phot[iph])
                    if j==0:
                        write_fits(s1nam0, ic.hdr, sp.Fnu_filt)
                    else:
                        write_fits(s1namj, ic.hdr, sp.Fnu_filt)
                    
        elif phot[iph]=='MIPS1' and Nch_LL>1:
            ##====
            ## LL
            ##====
            for i in trange(Nsub_LL,# leave=False,
                            desc=phot[iph]+' (LL synt phot) sub-maps'):
                f1nam = path_tmp+'LL/'+src+'_'+sub_LL[i]+'_LL'
                s1nam0 = cal_spec+'_'+sub_LL[i]+'_LL'

                for j in trange(Nmc+1-iresume, leave=False,
                                desc=' - '+sub_LL[i]+'_LL [MC]'):
                    j += iresume

                    s1namj = tmp_spec+'_'+sub_LL[i]+'_LL_'+str(j)
                    
                    ic = intercalib(f1nam+'_'+str(j))
                    sp = ic.synthetic_photometry(phot[iph])
                    if j==0:
                        write_fits(s1nam0, ic.hdr, sp.Fnu_filt)
                    else:
                        write_fits(s1namj, ic.hdr, sp.Fnu_filt)

    for i in range(Nsub_SL):
        s1nam0 = cal_spec+'_'+sub_SL[i]+'_SL'
        fits_Nmc = pathlib.Path(tmp_spec+'_'+sub_SL[i]+'_SL_'+str(Nmc)+fitsext)
        if Nmc>1 and fits_Nmc.exists():
            mcimage = []
            for j in range(Nmc):
                ds = read_fits(tmp_spec+'_'+sub_SL[i]+'_SL_'+str(j+1))
                mcimage.append(ds.data)
            mcimage = np.array(mcimage)
            unc = np.nanstd(mcimage, axis=0)
            write_fits(s1nam0+'_unc', ds.header, unc)

    for i in range(Nsub_LL):
        s1nam0 = cal_spec+'_'+sub_LL[i]+'_LL'
        fits_Nmc = pathlib.Path(tmp_spec+'_'+sub_LL[i]+'_LL_'+str(Nmc)+fitsext)
        if Nmc>1 and fits_Nmc.exists():
            mcimage = []
            for j in range(Nmc):
                ds = read_fits(tmp_spec+'_'+sub_LL[i]+'_LL_'+str(j+1))
                mcimage.append(ds.data)
            mcimage = np.array(mcimage)
            unc = np.nanstd(mcimage, axis=0)
            write_fits(s1nam0+'_unc', ds.header, unc)

## Sub-map intercalib correction (pixel/global scale)
##----------------------------------------------------
for iph in range(Nphot):
    print('\n Photometry: '+phot[iph])
    
    if photcal=='d' and phot[iph]=='IRAC4':
        cal_phot = path_cal+src+'_'+phot[iph]+'_DP'
    else:
        cal_phot = path_cal+src+'_'+phot[iph]+'_SINGS'
    cal_spec = path_cal+src+'_'+phot[iph]+'_IRS'
    tmp_spec = path_tmp+'calib/'+src+'_'+phot[iph]+'_IRS'
    
    if phot[iph]=='IRAC4' and Nch_SL>1:
        ##====
        ## SL
        ##====
        # xgrid = np.logspace(-2,4,10000)
        
        correct_sub = []
        correct_off = []
        correct_pixel = []
        for i in trange(Nsub_SL,# leave=False,
                        desc='Inter-calibrating IRS-SL sub-maps'):
            f1nam = path_tmp+'SL/'+src+'_'+sub_SL[i]+'_SL'
            p1nam0 = cal_phot+'_'+sub_SL[i]+'_SL'
            s1nam0 = cal_spec+'_'+sub_SL[i]+'_SL'

            for j in trange(Nmc+1, leave=False,
                            desc=' - '+src+'_'+sub_SL[i]+'_SL [MC]'):
                s1namj = tmp_spec+'_'+sub_SL[i]+'_SL_'+str(j)

                ## Data dictionary
                ## e.g. {subname0: val0, subname1: val1, ...}
                dict_spec = {}
                dict_spec_unc = {}
                if j==0:
                    dict_phot = {}
                    dict_phot_unc = {}
                    ds_phot = read_fits(p1nam0, p1nam0+'_unc')
                    Ny, Nx = ds_phot.data.shape
                    ds_spec = read_fits(s1nam0, s1nam0+'_unc')
                else:
                    ds_spec = read_fits(s1namj, s1nam0+'_unc') # intermediate unc

                for x in range(Nx):
                    for y in range(Ny):
                        subname = src+'_'+sub_SL[i]+'_SL_'+str(x+1)+'_'+str(y+1)
                        ## spec
                        dict_spec[subname] = ds_spec.data[y,x]
                        dict_spec_unc[subname] = ds_spec.unc[y,x]
                        if j==0:
                            ## Phot
                            dict_phot[subname] = ds_phot.data[y,x]
                            dict_phot_unc[subname] = ds_phot.unc[y,x]
                ## Dictionary to 1D array
                ## e.g. array([a1, a2, ...])
                pix_spec = np.array(list(dict_spec.values()))
                pix_spec_unc = np.array(list(dict_spec_unc.values()))
                if j==0:
                    pix_phot = np.array(list(dict_phot.values()))
                    pix_phot_unc = np.array(list(dict_phot_unc.values()))
                
                ## Mask NaNs and negtive values
                mask_nan = ~np.logical_or( np.isnan(pix_spec), np.isnan(pix_phot) )
                mask_neg = np.logical_and( pix_spec>0, pix_phot>0 )
                mask_lim = np.logical_and( pix_spec>1e-2, pix_phot>1e-2 )
                mask = np.logical_and( mask_nan, mask_neg)#, mask_lim )

                ## Add calibration error
                pix_phot_unc[i] = np.sqrt( pix_phot_unc[i]**2 + (pix_phot[i]*.03)**2 ) # IRAC 3% (Carey2010)
                pix_spec_unc[i]  = np.sqrt( pix_spec_unc[i]**2 + (pix_spec[i]*.05)**2 ) # IRS 5% (Carey2010)
                
                if mask.any():
                    
                    ## Linear fits (IRS phot - phot, SUB-MAPS)
                    ##=========================================
                    ## y = ax
                    popt0, pcov0 = curve_fit(f_lin0, pix_spec[mask], pix_phot[mask],
                                             sigma=pix_phot_unc[mask])
                    sub_gain0 = popt0[0]
                    ## y = ax + b
                    popt, pcov = curve_fit(f_lin, pix_spec[mask], pix_phot[mask],
                                           sigma=pix_phot_unc[mask])
                    sub_gain = popt[0]
                    sub_off = popt[1]
                    if j==0:
                        # print(src+'_'+sub_SL[i]+'_SL inter-calibration ('+phot[iph]
                        #       +') NO offset gain = {:.4}'.format(sub_gain0))

                        # print(src+'_'+sub_SL[i]+'_SL inter-calibration ('+phot[iph]
                        #       +') gain = {:.4}'.format(sub_gain))
                        # print(src+'_'+sub_SL[i]+'_SL inter-calibration ('+phot[iph]
                        #       +') offset = {:.4}'.format(sub_off))

                        xmin = min(pix_spec[mask])
                        xmax = max(pix_spec[mask])
                        ymin = min(pix_phot[mask])
                        ymax = max(pix_phot[mask])
                        lnxinf = np.log(xmin)
                        lnxsup = np.log(xmax)
                        # xgrid = np.logspace(lnxinf, lnxsup, 10000)
                        xgrid = np.linspace(xmin, ymax, 10000)

                        ## S/N ratio
                        # print(src+'_'+sub_SL[i]+'_SL S/N (IRS) = \n',
                        #       pix_spec[mask]/pix_spec_unc[mask])
                        # print(src+'_'+sub_SL[i]+'_SL S/N (phot) = \n',
                        #       pix_phot[mask]/pix_phot_unc[mask])

                        if sub_SL[i]=='04':
                            axlim = (1e1,1e4)
                            lab = '1'
                        elif sub_SL[i]=='06S':
                            axlim = (2e-3,1e4) # S
                            lab = '2'
                        elif sub_SL[i]=='06N':
                            axlim = (1e-2,2e3) # N
                            lab = '3'
                        elif sub_SL[i]=='08':
                            axlim = (5e-2,2e1)
                            lab = '4'
                        elif sub_SL[i]=='08c':
                            axlim = (None,None) # SL cap
                            lab = '-'
                        elif sub_SL[i]=='09N3':
                            axlim = (1e-1,2e1) # N3
                            lab = '5'
                        elif sub_SL[i]=='09N2':
                            axlim = (2e-1,5e1) # N2
                            lab = '6'

                        ## IRS phot - phot plot
                        p = pplot(pix_spec[mask], pix_phot[mask],
                                  yerr=pix_phot_unc[mask], xerr=pix_spec_unc[mask],
                                  fmt='s', ec='grey', elw=1, c='c',
                                  marker='.', markersize=20, capsize=2,
                                  xlog=1, ylog=1, nonposx='clip', nonposy='clip',
                                  xtkform='mylog', ytkform='log_sci',
                                  xlim=axlim, ylim=axlim,
                                  # title=phot[iph]+'_'+src+'_'+sub_SL[i]+'_SL',
                                  xlabel=r'$\rm IRS-sIRAC_{8\mu m}\ (MJy/sr)$',
                                  ylabel=r'$\rm IRAC_{8\mu m}\ (MJy/sr)$',
                                  figsize=(9,8),# right=.95, left=.15, bottom=.1, top=.95,
                                  loc='upper left',# anchor=(1,1),
                                  label=sub_SL[i]+'-SL',
                                  titlesize=20, xysize=20, tksize=20, legendsize=20)
                        p.add_plot(pix_spec[mask], pix_phot[mask],
                                   yerr=pix_phot_unc[mask], xerr=pix_spec_unc[mask],
                                   fmt='s', ec='grey', elw=.5, c='c', zorder=100,
                                   marker='.', markersize=.1, capsize=2,) # put errorbar on the front layer

                        label = 'y={0:.4}x'.format(sub_gain0)
                        p.add_plot(xgrid, f_lin0(xgrid, *popt0),
                                   c='k', ls='-', label=label, zorder=100,)

                        if sub_off<0:
                            label = 'y={0:.4}x{1:.4}'.format(sub_gain, sub_off)
                        else:
                            label = 'y={0:.4}x+{1:.4}'.format(sub_gain, sub_off)
                        p.add_plot(xgrid, f_lin(xgrid, *popt),
                                   c='k', ls='--', label=label, zorder=100,)
                        
                        # p.ax.text(.85,.05,'('+lab+')',size=30,c='grey',transform=p.ax.transAxes) # for the use of Hu_thesis
                        p.ax.legend(loc='upper left', fontsize=20, framealpha=0,)
                        
                        p.save(path_cal+'IC_'+phot[iph]+'_'+src+'_'+sub_SL[i]+'_SL.png',
                               transparent=True, figtight=True)
                
                ## Spectral correction
                ##=====================
                if j==0:
                    if write_SL=='y':
                    #     correct_sub.append(input(" - (SUB)MAP level correction (y/n)? "))
                    #     if correct_sub[i]=='y':
                    #         correct_off.append(input("   - Correct offset (y/n)? "))
                    #         correct_pixel.append(None)
                    #     else:
                    #         correct_off.append(None)
                    #         correct_pixel.append(input(" - PIXEL level correction (y/n)? "))
                    # else:
                    #     correct_sub.append(None)
                    #     correct_off.append(None)
                    #     correct_pixel.append(None)
                        correct_sub.append(None)
                        correct_off.append(None)
                        correct_pixel.append('y')
                
                if write_SL=='y':
                    pixel_gain = np.ones((Ny,Nx))
                    if mask.any():
                        for x in range(Nx):
                            for y in range(Ny):
                                subname = src+'_'+sub_SL[i]+'_SL_'+str(x+1)+'_'+str(y+1)
                                if ~ (np.isnan(dict_phot[subname]) or np.isnan(dict_spec[subname])):
                                    pixel_gain[y,x] *= dict_phot[subname]/dict_spec[subname]
                    if correct_pixel[i]=='y':
                        calib_gain = pixel_gain
                        calib_off = 0.
                    elif correct_sub[i]=='y':
                        if correct_off[i]=='y':
                            calib_gain = sub_gain
                            calib_off = sub_off
                        else:
                            calib_gain = sub_gain0
                            calib_off = 0.
                    else:
                        calib_gain = 1.0
                        calib_off = 0.
                    sc = intercalib(f1nam+'_'+str(j))
                    sc.correct_spec(calib_gain, calib_off, filOUT=f1nam+'_'+str(j))

    if phot[iph]=='IRAC4':
        for i in range(Nsub_SL):
            f1nam = path_tmp+'SL/'+src+'_'+sub_SL[i]+'_SL'
            if (write_SL=='y' and Nmc>1):
                fits_Nmc = pathlib.Path(f1nam+'_'+str(Nmc)+fitsext)
                if fits_Nmc.exists():
                    mcimage = []
                    for j in range(Nmc):
                        ds = read_fits(f1nam+'_'+str(j+1))
                        mcimage.append(ds.data)
                    mcimage = np.array(mcimage)
                    unc = np.nanstd(mcimage, axis=0)
                    write_fits(f1nam+'_unc', ds.header, unc, ds.wave)

        ##-----------------------
        ## Coadd all SL sub-maps
        ##-----------------------
        if coadd_SL=='y':
            if resume:
                iresume = int(input("Resume SL coadding from the iteration: "))
            else:
                iresume = 0
                
            for j in trange(Nmc+1-iresume,# leave=False,
                            desc='Coadding SL [MC]'):
                j += iresume
                
                fsub_SL = []
                for i in range(Nsub_SL):
                    f1nam = path_tmp+'SL/'+src+'_'+sub_SL[i]+'_SL'
                    fsub_SL.append(f1nam+'_'+str(j))
            
                if coadd_tool=='swarp':
                    swp = iswarp(refheader=coadd_footprint,
                                 tmpdir=path_tmp, verbose=verbose)
                    if j==0:
                        swp.combine(fsub_SL,# combtype='wgt_avg',
                                    keepedge=True,# cropedge=True,
                                    filOUT=path_out+src+'_SL')
                    else:
                        swp.combine(fsub_SL,# combtype='wgt_avg',
                                    keepedge=True,# cropedge=True,
                                    filOUT=path_tmp+src+'_SL_'+str(j))
                elif coadd_tool=='reproject':
                    mtg = imontage('exact', tmpdir=path_tmp)
                    if j==0:
                        mtg.coadd(fsub_SL, refheader=coadd_footprint,
                                  filOUT=path_out+src+'_SL')
                    else:
                        mtg.coadd(fsub_SL, refheader=coadd_footprint,
                                  filOUT=path_tmp+src+'_SL_'+str(j))

        ## Calculate unc
        fits_Nmc = pathlib.Path(path_tmp+src+'_SL_'+str(Nmc)+fitsext)
        if Nmc>1 and fits_Nmc.exists():
            mcimage = []
            for j in range(Nmc):
                ds = read_fits(path_tmp+src+'_SL_'+str(j+1))
                mcimage.append(ds.data)
            mcimage = np.array(mcimage)
            unc = np.nanstd(mcimage, axis=0)
            write_fits(path_out+src+'_SL_unc', ds.header, unc, ds.wave)

        ##------------------------
        ## Build SL slits (AKARI)
        ##------------------------
        if build_SL=='y':
            if resume:
                iresume = int(input("Resume SL slit building from the iteration: "))
            else:
                iresume = 0
                
            for iobs in trange(Nobs,# leave=False,
                              desc='Reproject to IRC slits'):
                f2nam = path_tmp+'SL/'+parobs[iobs][0]
                ## refheader
                refheader = fixwcs(out_irc[iobs]+fitsext).header
                
                for j in trange(Nmc+1-iresume, leave=False,
                                desc='Coadding SL [MC]'):
                    j += iresume
                    
                    fsub_SL = []
                    for i in range(Nsub_SL):
                        f1nam = path_tmp+'SL/'+src+'_'+sub_SL[i]+'_SL'
                        fsub_SL.append(f1nam+'_'+str(j))
                
                    if coadd_tool=='swarp':
                        swp = iswarp(refheader=refheader,
                                     tmpdir=path_tmp+'SL/', verbose=verbose)
                        swp.combine(fsub_SL,# combtype='wgt_avg',
                                    keepedge=True,# cropedge=True,
                                    filOUT=f2nam+'_MC_'+str(j))
                    elif coadd_tool=='reproject':
                        mtg = imontage('exact', tmpdir=path_tmp)
                        mtg.coadd(fsub_SL, refheader=refheader,
                                  filOUT=f2nam+'_MC_'+str(j))

        for iobs in range(Nobs):
            f2nam = path_tmp+'SL/'+parobs[iobs][0]
            ## Calculate unc
            fits_Nmc = pathlib.Path(f2nam+'_MC_'+str(Nmc)+fitsext)
            if Nmc>1 and fits_Nmc.exists():
                mcimage = []
                for j in range(Nmc):
                    ds = read_fits(f2nam+'_MC_'+str(j+1))
                    mcimage.append(ds.data)
                mcimage = np.array(mcimage)
                unc = np.nanstd(mcimage, axis=0)
                write_fits(f2nam+'_MC_unc', ds.header, unc, ds.wave)

    if phot[iph]=='MIPS1' and Nch_LL>1:
        ##====
        ## LL
        ##====
        # xgrid = np.logspace(-2,4,10000)
        
        correct_sub = []
        correct_off = []
        correct_pixel = []
        for i in trange(Nsub_LL,# leave=False,
                        desc='Inter-calibrating IRS-LL sub-maps'):
            f1nam = path_tmp+'LL/'+src+'_'+sub_LL[i]+'_LL'
            p1nam0 = cal_phot+'_'+sub_LL[i]+'_LL'
            s1nam0 = cal_spec+'_'+sub_LL[i]+'_LL'

            for j in trange(Nmc+1, leave=False,
                            desc=' - '+src+'_'+sub_LL[i]+'_LL [MC]'):
                s1namj = tmp_spec+'_'+sub_LL[i]+'_LL_'+str(j)
                
                ds_phot = read_fits(p1nam0, p1nam0+'_unc')
                Ny, Nx = ds_phot.data.shape
                if j==0:
                    ds_spec = read_fits(s1nam0, s1nam0+'_unc')
                else:
                    ds_spec = read_fits(s1namj, s1nam0+'_unc') # intermediate unc
                
                ## Data dictionary
                ## e.g. {subname0: val0, subname1: val1, ...}
                dict_phot = {}
                dict_phot_unc = {}
                dict_spec = {}
                dict_spec_unc = {}
                for x in range(Nx):
                    for y in range(Ny):
                        subname = src+'_'+sub_LL[i]+'_LL_'+str(x+1)+'_'+str(y+1)
                        ## Phot
                        dict_phot[subname] = ds_phot.data[y,x]
                        dict_phot_unc[subname] = ds_phot.unc[y,x]
                        ## spec
                        dict_spec[subname] = ds_spec.data[y,x]
                        dict_spec_unc[subname] = ds_spec.unc[y,x]
                
                ## Dictionary to 1D array
                ## e.g. array([a1, a2, ...])
                pix_phot = np.array(list(dict_phot.values()))
                pix_phot_unc = np.array(list(dict_phot_unc.values()))
                pix_spec = np.array(list(dict_spec.values()))
                pix_spec_unc = np.array(list(dict_spec_unc.values()))
                
                ## Mask NaNs and negtive values
                mask_nan = ~np.logical_or( np.isnan(pix_spec), np.isnan(pix_phot) )
                mask_neg = np.logical_and( pix_spec>0, pix_phot>0 )
                mask_lim = np.logical_and( pix_spec>1e-2, pix_phot>1e-2 )
                mask = np.logical_and( mask_nan, mask_neg,)# mask_lim )

                ## Add calibration error
                pix_phot_unc[i] = np.sqrt( pix_phot_unc[i]**2 + (pix_phot[i]*.04)**2 ) # MIPS 4% (Carey2010)
                pix_spec_unc[i]  = np.sqrt( pix_spec_unc[i]**2 + (pix_spec[i]*.05)**2 ) # IRS 5% (Carey2010)
                
                if mask.any():
                    
                    ## Linear fits (IRS phot - phot, SUB-MAPS)
                    ## For M82, only SINGS has MIPS1
                    ##=========================================
                    ## y = ax
                    popt0, pcov0 = curve_fit(f_lin0, pix_spec[mask], pix_phot[mask],
                                             sigma=pix_phot_unc[mask])
                    sub_gain0 = popt0[0]
                    ## y = ax + b
                    popt, pcov = curve_fit(f_lin, pix_spec[mask], pix_phot[mask],
                                           sigma=pix_phot_unc[mask])
                    sub_gain = popt[0]
                    sub_off = popt[1]
                    if j==0:
                        # print(src+'_'+sub_LL[i]+'_LL inter-calibration ('+phot[iph]
                        #       +') NO offset gain = {:.4}'.format(sub_gain0))

                        # print(src+'_'+sub_LL[i]+'_LL inter-calibration ('+phot[iph]
                        #       +') gain = {:.4}'.format(sub_gain))
                        # print(src+'_'+sub_LL[i]+'_LL inter-calibration ('+phot[iph]
                        #       +') offset = {:.4}'.format(sub_off))

                        xmin = min(pix_spec[mask])
                        xmax = max(pix_spec[mask])
                        ymin = min(pix_phot[mask])
                        ymax = max(pix_phot[mask])
                        lnxinf = np.log(xmin)
                        lnxsup = np.log(xmax)
                        # xgrid = np.logspace(lnxinf, lnxsup, 10000)
                        xgrid = np.linspace(xmin, ymax, 10000)

                        ## S/N ratio
                        # print(src+'_'+sub_LL[i]+'_LL S/N (IRS) = \n',
                        #       pix_spec[mask]/pix_spec_unc[mask])
                        # print(src+'_'+sub_LL[i]+'_LL S/N (SINGS) = \n',
                        #       pix_phot[mask]/pix_phot_unc[mask])

                        if sub_LL[i]=='04':
                            axlim = (1e0,1e4)
                            lab = 'a'
                        elif sub_LL[i]=='05':
                            axlim = (8e-1,5e2)
                            lab = '7'
                        elif sub_LL[i]=='06':
                            axlim = (2e-4,2e4)
                            lab = '8'
                        elif sub_LL[i]=='08':
                            axlim = (2e-5,5e1)
                            lab = '9'
                        elif sub_LL[i]=='09N3':
                            axlim = (5e-5,1e3) # N3
                            lab = '10'
                        elif sub_LL[i]=='09N5':
                            axlim = (8e-2,1e1) # N5
                            lab = '11'
                        elif sub_LL[i]=='09N2':
                            axlim = (1e-1,2e1) # N2
                            lab = '12'
                    
                        ## IRS - phot plot (For M82, only SINGS has MIPS1)
                        p = pplot(pix_spec[mask], pix_phot[mask],
                                  yerr=pix_phot_unc[mask], xerr=pix_spec_unc[mask],
                                  fmt='s', ec='grey', elw=1, c='c',
                                  marker='.', markersize=20, capsize=2,
                                  xlog=1, ylog=1, nonposx='clip', nonposy='clip',
                                  xtkform='mylog', ytkform='log_sci',
                                  xlim=axlim, ylim=axlim,
                                  # title=phot[iph]+'_'+src+'_'+sub_LL[i]+'_LL',
                                  xlabel=r'$\rm IRS-sMIPS_{24\mu m}\ (MJy/sr)$',
                                  ylabel=r'$\rm MIPS_{24\mu m}\ (MJy/sr)$',
                                  # figsize=(11,8), right=.78, left=.12, bottom=.1, top=.95,
                                  figsize=(9,8),# right=.95, left=.15, bottom=.1, top=.95,
                                  loc='upper left',# anchor=(1,1),
                                  label=sub_LL[i]+'-LL',
                                  titlesize=20, xysize=20, tksize=20, legendsize=20)
                        p.add_plot(pix_spec[mask], pix_phot[mask],
                                   yerr=pix_phot_unc[mask], xerr=pix_spec_unc[mask],
                                   fmt='s', ec='grey', elw=.5, c='c', zorder=100,
                                   marker='.', markersize=.1, capsize=2,) # put errorbar on the front layer

                        label = 'y={0:.4}x'.format(sub_gain0)
                        p.add_plot(xgrid, f_lin0(xgrid, *popt0),
                                   c='k', ls='-', label=label, zorder=100,)

                        if sub_off<0:
                            label = 'y={0:.4}x{1:.4}'.format(sub_gain, sub_off)
                        else:
                            label = 'y={0:.4}x+{1:.4}'.format(sub_gain, sub_off)
                        p.add_plot(xgrid, f_lin(xgrid, *popt),
                                   c='k', ls='--', label=label, zorder=100,)
                        
                        # p.ax.text(.85,.05,'('+lab+')',size=30,c='grey',transform=p.ax.transAxes) # for the use of Hu_thesis
                        p.ax.legend(loc='upper left', fontsize=20, framealpha=0,)
                        
                        p.save(path_cal+'IC_'+phot[iph]+'_'+src+'_'+sub_LL[i]+'_LL.png',
                               transparent=True, figtight=True)
                    
                ## Spectral correction
                ##=====================
                if j==0:
                    if write_LL=='y':
                    #     correct_sub.append(input(" - (SUB)MAP level correction (y/n)? "))
                    #     if correct_sub[i]=='y':
                    #         correct_off.append(input("   - Correct offset (y/n)? "))
                    #         correct_pixel.append(None)
                    #     else:
                    #         correct_off.append(None)
                    #         correct_pixel.append(input(" - PIXEL level correction (y/n)? "))
                    # else:
                    #     correct_sub.append(None)
                    #     correct_off.append(None)
                    #     correct_pixel.append(None)
                        correct_sub.append(None)
                        correct_off.append(None)
                        correct_pixel.append('y')
                
                if write_LL=='y':
                    pixel_gain = np.ones((Ny,Nx))
                    if mask.any():
                        for x in range(Nx):
                            for y in range(Ny):
                                subname = src+'_'+sub_LL[i]+'_LL_'+str(x+1)+'_'+str(y+1)
                                if ~ (np.isnan(dict_phot[subname]) or np.isnan(dict_spec[subname])):
                                    pixel_gain[y,x] *= dict_phot[subname]/dict_spec[subname]
                    if correct_pixel[i]=='y':
                        calib_gain = pixel_gain
                        calib_off = 0.
                    elif correct_sub[i]=='y':
                        if correct_off[i]=='y':
                            calib_gain = sub_gain
                            calib_off = sub_off
                        else:
                            calib_gain = sub_gain0
                            calib_off = 0.
                    else:
                        calib_gain = 1.0
                        calib_off = 0.
                    sc = intercalib(f1nam+'_'+str(j))
                    sc.correct_spec(calib_gain, calib_off, filOUT=f1nam+'_'+str(j))

    if phot[iph]=='MIPS1':
        for i in range(Nsub_LL):
            f1nam = path_tmp+'LL/'+src+'_'+sub_LL[i]+'_LL'
            if (write_LL=='y' and Nmc>1):
                fits_Nmc = pathlib.Path(f1nam+'_'+str(Nmc)+fitsext)
                if fits_Nmc.exists():
                    mcimage = []
                    for j in range(Nmc):
                        ds = read_fits(f1nam+'_'+str(j+1))
                        mcimage.append(ds.data)
                    mcimage = np.array(mcimage)
                    unc = np.nanstd(mcimage, axis=0)
                    write_fits(f1nam+'_unc', ds.header, unc, ds.wave)

        ##-----------------------
        ## Coadd all LL sub-maps
        ##-----------------------
        if coadd_LL=='y':
            if resume:
                iresume = int(input("Resume LL coadding from the iteration: "))
            else:
                iresume = 0
                
            for j in trange(Nmc+1-iresume, leave=False,
                            desc='Coadding LL [MC]'):
                j += iresume
                
                fsub_LL = []
                for i in range(Nsub_LL):
                    f1nam = path_tmp+'LL/'+src+'_'+sub_LL[i]+'_LL'
                    fsub_LL.append(f1nam+'_'+str(j))
            
                if coadd_tool=='swarp':
                    swp = iswarp(refheader=coadd_footprint,
                                 tmpdir=path_tmp, verbose=verbose)
                    if j==0:
                        swp.combine(fsub_LL,# combtype='wgt_avg',
                                    keepedge=True,# cropedge=True,
                                    filOUT=path_out+src+'_LL')
                    else:
                        swp.combine(fsub_LL,# combtype='wgt_avg',
                                    keepedge=True,# cropedge=True,
                                    filOUT=path_tmp+src+'_LL_'+str(j))
                elif coadd_tool=='reproject':
                    mtg = imontage('exact', tmpdir=path_tmp)
                    if j==0:
                        mtg.coadd(fsub_LL, refheader=coadd_footprint,
                                  filOUT=path_out+src+'_LL')
                    else:
                        mtg.coadd(fsub_LL, refheader=coadd_footprint,
                                  filOUT=path_tmp+src+'_LL_'+str(j))

            ## Calculate unc
            fits_Nmc = pathlib.Path(path_tmp+src+'_LL_'+str(Nmc)+fitsext)
            if Nmc>1 and fits_Nmc.exists():
                mcimage = []
                for j in range(Nmc):
                    ds = read_fits(path_tmp+src+'_LL_'+str(j+1))
                    mcimage.append(ds.data)
                mcimage = np.array(mcimage)
                unc = np.nanstd(mcimage, axis=0)
                write_fits(path_out+src+'_LL_unc', ds.header, unc, ds.wave)

        ##------------------------
        ## Build LL slits (AKARI)
        ##------------------------
        if build_LL=='y':
            if resume:
                iresume = int(input("Resume LL slit building from the iteration: "))
            else:
                iresume = 0
                
            for iobs in trange(Nobs,# leave=False,
                              desc='Reproject to IRC slits'):
                f2nam = path_tmp+'LL/'+parobs[iobs][0]
                ## refheader
                refheader = fixwcs(out_irc[iobs]+fitsext).header
                
                for j in trange(Nmc+1-iresume, leave=False,
                                desc='Coadding LL [MC]'):
                    j += iresume
                    
                    fsub_LL = []
                    for i in range(Nsub_LL):
                        f1nam = path_tmp+'LL/'+src+'_'+sub_LL[i]+'_LL'
                        fsub_LL.append(f1nam+'_'+str(j))
                
                    if coadd_tool=='swarp':
                        swp = iswarp(refheader=refheader,
                                     tmpdir=path_tmp+'LL/', verbose=verbose)
                        swp.combine(fsub_LL,# combtype='wgt_avg',
                                    keepedge=True,# cropedge=True,
                                    filOUT=f2nam+'_MC_'+str(j))
                    elif coadd_tool=='reproject':
                        mtg = imontage('exact', tmpdir=path_tmp)
                        mtg.coadd(fsub_LL, refheader=refheader,
                                  filOUT=f2nam+'_MC_'+str(j))

            for iobs in range(Nobs):
                f2nam = path_tmp+'LL/'+parobs[iobs][0]
                ## Calculate unc
                fits_Nmc = pathlib.Path(f2nam+'_MC_'+str(Nmc)+fitsext)
                if Nmc>1 and fits_Nmc.exists():
                    mcimage = []
                    for j in range(Nmc):
                        ds = read_fits(f2nam+'_MC_'+str(j+1))
                        mcimage.append(ds.data)
                    mcimage = np.array(mcimage)
                    unc = np.nanstd(mcimage, axis=0)
                    write_fits(f2nam+'_MC_unc', ds.header, unc, ds.wave)


print('\n>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>')
print('\nSpitzer IRS => AKARI IRC\n')
print('>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n')


##----------------------------------------------------------

##              Inter-calibration (global)

##----------------------------------------------------------
pre_calib0 = input("Prepare photometry (y/n)? ")
synt_phot0 = input("Prepare IRS synthetic photometry (y/n)? ")
do_intcal0 = input("Do inter-calibration (y/n)? ")
if (synt_phot0!='y' or pre_calib0!='y'):
    warnings.warn('You need to prepare (synthetic) photometry if you have not done it yet.')

os.makedirs(path_tmp+'calib/', exist_ok=True)

## Homogenize photometry and synthetic photometry
##------------------------------------------------
for iph in range(Nphot):
    
    ## Prepare photometry
    if phot[iph]=='IRAC4':
        phot_ker = path_ker+'Kernel_HiRes_IRAC_8.0_to_Gauss_06.0'
    elif phot[iph]=='MIPS1':
        phot_ker = path_ker+'Kernel_HiRes_MIPS_24_to_Gauss_06.0'

    if photcal=='d' and phot[iph]=='IRAC4':
        cal_phot = path_cal+src+'_'+phot[iph]+'_DP'
        tmp_phot = path_tmp+'calib/'+src+'_'+phot[iph]+'_DP'
    else:
        cal_phot = path_cal+src+'_'+phot[iph]+'_SINGS'
        tmp_phot = path_tmp+'calib/'+src+'_'+phot[iph]+'_SINGS'
    cal_spec = path_cal+src+'_'+phot[iph]+'_IRS'
    tmp_spec = path_tmp+'calib/'+src+'_'+phot[iph]+'_IRS'
    
    ##--------------
    ## SINGS (phot)
    ##--------------
    if pre_calib0=='y':
        if photcal=='d' and phot[iph]=='IRAC4':
            raw_phot = path_phot+src+'_'+phot[iph]+'_DP'
        else:
            raw_phot = path_phot+src+'_'+phot[iph]+'_SINGS'

        ## Reproject phot (atlas)
        if coadd_tool=='swarp':
            swp = iswarp(refheader=coadd_footprint, tmpdir=path_tmp)
            swp.combine_mc(raw_phot, keepedge=True, #cropedge=True,
                           dist='norm', Nmc=Nmc, filOUT=tmp_phot+'_MC')
        elif coadd_tool=='reproject':
            mtg = imontage('exact', tmpdir=path_tmp)
            mtg.reproject_mc(raw_phot, refheader=coadd_footprint,
                             dist='norm', Nmc=Nmc, filOUT=tmp_phot+'_MC')

    for iobs in trange(Nobs, #leave=False,
                       desc=phot[iph]+' processing'):
        xscale, yscale = read_hdf5(filog+parobs[iobs][0], 'Super pixel size')
        if phot[iph]=='IRAC4':
            f2nam = path_tmp+'SL/'+parobs[iobs][0]
        elif phot[iph]=='MIPS1':
            f2nam = path_tmp+'LL/'+parobs[iobs][0]
        p2nam0 = cal_phot+'_'+parobs[iobs][0]
        s2nam0 = cal_spec+'_'+parobs[iobs][0]

        fits_Nmc0 = pathlib.Path(f2nam+'_MC_0'+fitsext)
        if pre_calib0=='y' and fits_Nmc0.exists():
            refheader = fixwcs(f2nam+'_MC_0'+fitsext).header
            
            for j in trange(Nmc+1, leave=False,
                            desc=' - '+parobs[iobs][0]+' [MC]'):
                p2namj = tmp_phot+'_'+parobs[iobs][0]+'_'+str(j)
                
                if iobs==0:
                    ## Convolve phot
                    if j==0:
                        conv = iconvolve(tmp_phot+'_MC',
                                         kfile=phot_ker, klist=csv_ker,
                                         filOUT=tmp_phot+'_MC')
                    else:
                        conv = iconvolve(tmp_phot+'_MC_'+str(j),
                                         kfile=phot_ker, klist=csv_ker,
                                         filOUT=tmp_phot+'_MC_'+str(j))
                    conv.do_conv(idldir=path_idl)
        
                ## Reproject phot (slit)
                if coadd_tool=='swarp':
                    swp = iswarp(refheader=refheader, tmpdir=path_tmp)
                    if j==0:
                        swp.combine(tmp_phot+'_MC', filOUT=p2nam0)
                    else:
                        swp.combine(tmp_phot+'_MC_'+str(j), filOUT=p2namj)
                elif coadd_tool=='reproject':
                    mtg = imontage('exact', tmpdir=path_tmp)
                    if j==0:
                        mtg.reproject(tmp_phot+'_MC', refheader=refheader,
                                      filOUT=p2nam0)
                    else:
                        mtg.reproject(tmp_phot+'_MC_'+str(j), refheader=refheader,
                                      filOUT=p2namj)
                
                if j==0:
                    igroupixel(p2nam0,
                               xscale=xscale, yscale=yscale,
                               filOUT=p2nam0)
                else:
                    igroupixel(p2namj,
                               xscale=xscale, yscale=yscale,
                               filOUT=p2namj)
            ## Calculate unc
            fits_Nmc = pathlib.Path(tmp_phot+'_'+parobs[iobs][0]+'_'+str(Nmc)+fitsext)
            if Nmc>1 and fits_Nmc.exists():
                mcimage = []
                for j in range(Nmc):
                    ds = read_fits(tmp_phot+'_'+parobs[iobs][0]+'_'+str(j+1))
                    mcimage.append(ds.data)
                mcimage = np.array(mcimage)
                unc = np.nanstd(mcimage, axis=0)
                write_fits(p2nam0+'_unc', ds.header, unc)

        ##-----------------------------
        ## Synthetic photometry (spec)
        ##-----------------------------
        if synt_phot0=='y':# and np.logical_or(phot[iph]=='IRAC4' and Nch_SL>1,
                                             # phot[iph]=='MIPS1' and Nch_LL>1):
            for j in trange(Nmc+1, leave=False,
                            desc=' - '+parobs[iobs][0]+' (IRS sythetic photometry) [MC]'):
                s2namj = tmp_spec+'_'+parobs[iobs][0]+'_'+str(j)
                
                ic = intercalib(f2nam+'_MC_'+str(j))
                sp = ic.synthetic_photometry(phot[iph], xscale=xscale, yscale=yscale)
                if j==0:
                    write_fits(s2nam0, ic.hdr, sp.Fnu_filt)
                else:
                    write_fits(s2namj, ic.hdr, sp.Fnu_filt)
        
            ## Calculate unc
            fits_Nmc = pathlib.Path(tmp_spec+'_'+parobs[iobs][0]+'_'+str(Nmc)+fitsext)
            if Nmc>1 and fits_Nmc.exists():
                mcimage = []
                for j in range(Nmc):
                    ds = read_fits(tmp_spec+'_'+parobs[iobs][0]+'_'+str(j+1))
                    mcimage.append(ds.data)
                mcimage = np.array(mcimage)
                unc = np.nanstd(mcimage, axis=0)
                write_fits(cal_spec+'_'+parobs[iobs][0]+'_unc', ds.header, unc)

## Intercalib correction (pixel/global scale)
##--------------------------------------------
for iph in range(Nphot):
    if phot[iph]=='IRAC4':
        mod = 'SL'
    elif phot[iph]=='MIPS1':
        mod = 'LL'
    else:
        mod = None

    write_irs = input("Write IRS "+mod+" spectra (y/n)? ")

    # if phot[iph]=='IRAC4' and Nch_SL>1:
    #     do_correct = True
    # elif phot[iph]=='MIPS1' and Nch_LL>1:
    #     do_correct = True
    # else:
    #     do_correct = False

    # if do_correct:
    print('\n Photometry: '+phot[iph])
    
    if photcal=='d' and phot[iph]=='IRAC4':
        cal_phot = path_cal+src+'_'+phot[iph]+'_DP'
    else:
        cal_phot = path_cal+src+'_'+phot[iph]+'_SINGS'
    cal_spec = path_cal+src+'_'+phot[iph]+'_IRS'
    tmp_spec = path_tmp+'calib/'+src+'_'+phot[iph]+'_IRS'
    
    for j in trange(Nmc+1, #leave=False,
                    desc='IRS spectral correction [MC]'):
        ## Data dictionary
        ## e.g. [{'A0': a0, 'A1': a1, ...}, {'B0': b0, 'B1': b1, ...}, ...]
        dict_spec = []
        dict_spec_unc = []
        ## Data 1D array
        ## e.g. [array([a1, a2, ...]), array([b1, b2, ...]), ...]
        pix_spec = []
        pix_spec_unc = []
        if j==0:
            dict_phot = []
            dict_phot_unc = []
            pix_phot = []
            pix_phot_unc = []
            pix_mask_fit = []
            # xgrid = np.logspace(-2,4,10000)
    
            xmin = []
            xmax = []
            ymin = []
            ymax = []
        for iobs in range(Nobs):
            if phot[iph]=='IRAC4':
                f2nam = path_tmp+'SL/'+parobs[iobs][0]
            elif phot[iph]=='MIPS1':
                f2nam = path_tmp+'LL/'+parobs[iobs][0]
            p2nam0 = cal_phot+'_'+parobs[iobs][0]
            s2nam0 = cal_spec+'_'+parobs[iobs][0]
            s2namj = tmp_spec+'_'+parobs[iobs][0]+'_'+str(j)
            
            # Nx = read_hdf5(filog+parobs[iobs][0], 'Slit width')
            Ny = read_hdf5(filog+parobs[iobs][0], 'Slit length')[0]
            xscale, yscale = read_hdf5(filog+parobs[iobs][0], 'Super pixel size')
            # Nxs = math.ceil(Nx/xscale)
            Nys = math.ceil(Ny/yscale)
    
            data_spec = {}
            unc_spec = {}
            mask0 = []
            if j==0:
                data_phot = {}
                unc_phot = {}
                ds_phot = read_fits(p2nam0, p2nam0+'_unc')
                ds_spec = read_fits(s2nam0, s2nam0+'_unc')
            else:
                ds_spec = read_fits(s2namj, s2nam0+'_unc') # intermediate unc
            
            for ys in range(Nys):
                subname = parobs[iobs][0][0]+str(ys+1)
                yarr = sup2pix(ys, yscale, Npix=Ny, origin=0)
                ## spec
                data_spec[subname] = ds_spec.data[yarr[0],0]
                unc_spec[subname] = ds_spec.unc[yarr[0],0]
                if j==0:
                    ## Phot
                    data_phot[subname] = ds_phot.data[yarr[0],0]
                    unc_phot[subname] = ds_phot.unc[yarr[0],0]
                    
                ## Do not fit pixels affected by the MIPS1 map hole
                if phot[iph]=='MIPS1' and \
                    (subname=='A3' or subname=='A4' or subname=='A5' or subname=='A6' \
                     or subname=='A7' or subname=='B4' or subname=='B5'):
                    mask0.append(False)
                else:
                    mask0.append(True)
            dict_spec.append(data_spec)
            dict_spec_unc.append(unc_spec)
            ## Dictionary to 1D array
            pix_spec.append(np.array(list(dict_spec[iobs].values())))
            pix_spec_unc.append(np.array(list(dict_spec_unc[iobs].values())))
            if j==0:
                dict_phot.append(data_phot)
                dict_phot_unc.append(unc_phot)
                pix_phot.append(np.array(list(dict_phot[iobs].values())))
                pix_phot_unc.append(np.array(list(dict_phot_unc[iobs].values())))
                pix_mask_fit.append(mask0)
        
            ## Mask NaNs and negtive values
            mask_nan = ~np.logical_or( np.isnan(pix_spec[iobs]), np.isnan(pix_phot[iobs]) )
            mask_neg = np.logical_and( pix_spec[iobs]>0, pix_phot[iobs]>0 )
            mask_lim = np.logical_and( pix_spec[iobs]>1e-2, pix_phot[iobs]>1e-2 )
            mask = np.logical_and( mask_nan, mask_neg)#, mask_lim )
            maskfit = np.logical_and( mask, mask0)#, mask_lim )
    
            ## Add calibration error
            if phot[iph]=='IRAC4' and Nch_SL>1:
                ylabel = r'$\rm IRAC_{8\mu m}\ (MJy/sr)$'
                # lab = '(c)' # for the use of Hu_thesis
                pix_phot_unc[iobs]  = np.sqrt( pix_phot_unc[iobs]**2 + (pix_phot[iobs]*.03)**2 ) # IRAC 3% (Carey2010)
            elif phot[iph]=='MIPS1' and Nch_LL>1:
                ylabel = r'$\rm MIPS_{24\mu m}\ (MJy/sr)$'
                # lab = '(d)' # for the use of Hu_thesis
                pix_phot_unc[iobs] = np.sqrt( pix_phot_unc[iobs]**2 + (pix_phot[iobs]*.04)**2 ) # MIPS 4% (Carey2010)
            pix_spec_unc[iobs] = np.sqrt( pix_spec_unc[iobs]**2 + (pix_spec[iobs]*.05)**2 ) # IRS 5% (Carey2010)
    
            if mask.any():
                
                if j==0:
                    ## Linear fit (IRS phot - phot, SLITS)
                    ##=====================================
                    popt, pcov = curve_fit(f_lin0, pix_spec[iobs][maskfit], pix_phot[iobs][maskfit],
                                           sigma=pix_phot_unc[iobs][maskfit])
                    slit_gain = popt[0]
                    # slit_gain = 1.
                    # slit_off = popt[1]
                    slit_off = 0.
                    # slit_off = popt[0]
                    # print(parobs[iobs][0]+' inter-Calibration ('+phot[iph]+
                    #       ') offset = {:.4}'.format(slit_off))
                        
                    ## S/N ratio
                    # print(parobs[iobs][0]+' S/N (IRS) = \n',
                    #       pix_spec[iobs][mask]/pix_spec_unc[iobs][mask])
                    # print(parobs[iobs][0]+' S/N (phot) = \n',
                    #       pix_phot[iobs][mask]/pix_phot_unc[iobs][mask])
                        
                    if do_intcal0=='y':
                        
                        ## IRS - phot plot
                        if iobs==0:
                            p = pplot(fmt='s', ec='grey', elw=1,
                                      xlog=1, ylog=1, nonposx='clip', nonposy='clip',
                                      xtkform='mylog', ytkform='log_sci',
                                      # xlim=(1e-2,1e4), ylim=(1e-2,1e4),
                                      # title=src+'_'+phot[iph],
                                      xlabel='IRS (MJy/sr)', ylabel=ylabel,
                                      figsize=(11,8), right=.78, left=.12, bottom=.1, top=.95,
                                      loc='upper left', anchor=(1,1),
                                      titlesize=20, xysize=20, tksize=20, legendsize=20)
                            
                        p.add_plot(pix_spec[iobs][mask], pix_phot[iobs][mask],
                                   yerr=pix_phot_unc[iobs][mask], xerr=pix_spec_unc[iobs][mask],
                                   fmt='s', ec='grey', c=colors[iobs+1],
                                   marker=markers[iobs], markersize=20, capsize=2,
                                   label=parobs[iobs][0])
                        p.add_plot(pix_spec[iobs][mask], pix_phot[iobs][mask],
                                   yerr=pix_phot_unc[iobs][mask], xerr=pix_spec_unc[iobs][mask],
                                   fmt='s', ec='grey', c=colors[iobs+1], zorder=100,
                                   marker=markers[iobs], markersize=.1, capsize=2) # put errorbar on the front layer
                        
                        p.save(path_cal+'IC_'+phot[iph]+'.png', transparent=True, figtight=True)
    
                        xmin.append(min(pix_spec[iobs][mask]))
                        xmax.append(max(pix_spec[iobs][mask]))
                        ymin.append(min(pix_phot[iobs][mask]))
                        ymax.append(max(pix_phot[iobs][mask]))
                        lnxinf = np.log(min(xmin))
                        lnxsup = np.log(max(xmax))
                        # xgrid = np.logspace(lnxinf, lnxsup, 10000)
                        xgrid = np.linspace(xmin, xmax, 100000)
                        
        ## Linear fit (IRS phot - phot, ATLAS)
        ##=====================================
        fitx = np.concatenate(pix_spec, axis=0)
        uncx = np.concatenate(pix_spec_unc, axis=0)
        fity = np.concatenate(pix_phot, axis=0)
        uncy = np.concatenate(pix_phot_unc, axis=0)
        mask0 = np.concatenate(pix_mask_fit, axis=0)
        
        mask_nan = ~np.logical_or( np.isnan(fitx), np.isnan(fity) )
        mask_neg = np.logical_and( fitx>0, fity>0 )
        mask_lim = np.logical_and( fitx>1e-2, fity>1e-2 )
        mask = np.logical_and( mask_nan, mask_neg )
        maskfit = np.logical_and( mask, mask0)#, mask_lim )
        
        popt, pcov = curve_fit(f_lin0, fitx[maskfit], fity[maskfit],
                               sigma=uncy[maskfit])
        atlas_gain = popt[0]
        # atlas_gain = 1.
        # atlas_off = popt[1]
        atlas_off = 0.
        # atlas_off = popt[0]
        
        if j==0:
            if do_intcal0=='y':
                print('Atlas inter-calibration ('+phot[iph]+') gain = {:.4}'.format(atlas_gain))
                # print('Atlas inter-Calibration ('+phot[iph]+') offset = {:.4}'.format(atlas_off))
                
                label = 'y={0:.4}x'.format(atlas_gain)
                if iph==1:
                    label = 'y={0:.2}x'.format(atlas_gain)
                # label = 'y={0:.4}x+{1:.4}'.format(atlas_gain, atlas_off)
                p.add_plot(xgrid, f_lin0(xgrid, *popt),
                           c='k', ls='-', label=label)
                # p.ax.text(.9,.05,lab,size=30,c='grey',
                #           transform=p.ax.transAxes) # for the use of Hu_thesis
                p.ax.legend(loc='upper left', bbox_to_anchor=(1,1),
                            fontsize=20, framealpha=0,)
                
                p.save(path_cal+'IC_'+phot[iph]+'.png',
                       transparent=True, figtight=True)
    
        ## Spectral correction
        ##=====================
        if j==0:
            if write_irs=='y':
                # if np.logical_or(phot[iph]=='IRAC4' and build_SL=='y',
                #                  phot[iph]=='MIPS' and build_LL=='y'):
                correct_atlas = input(" - ATLAS level correction (y/n)? ")
                if correct_atlas=='y':
                    correct_pixel = None
                else:
                    correct_pixel = input(" - PIXEL level correction (y/n)? ")
                # else:
                #     warnings.warn('Rebuild slits before spectral correction.')
                #     correct_atlas = None
                #     correct_pixel = None
            else:
                correct_atlas = None
                correct_pixel = None
        
        for iobs in range(Nobs):
            if phot[iph]=='IRAC4':
                f2nam = path_tmp+'SL/'+parobs[iobs][0]
            elif phot[iph]=='MIPS1':
                f2nam = path_tmp+'LL/'+parobs[iobs][0]
            xscale, yscale = read_hdf5(filog+parobs[iobs][0], 'Super pixel size')
            mask_nan = ~np.logical_or( np.isnan(pix_spec[iobs]), np.isnan(pix_phot[iobs]) )
            mask_neg = np.logical_and( pix_spec[iobs]>0, pix_phot[iobs]>0 )
            mask = np.logical_and( mask_nan, mask_neg )
            sc = intercalib(f2nam+'_MC_'+str(j))
            Nw, Ny, Nx = sc.im.shape
            pixel_gain = np.ones((Ny,Nx))
            if mask.any():
                for y in range(Ny):
                    ys = pix2sup(y, yscale, origin=0)
                    if mask[ys]:
                        pixel_gain[y,:] *= pix_phot[iobs][ys]/pix_spec[iobs][ys]
                    # if y%yscale==0:
                    #     print(parobs[iobs][0]+' inter-calibration gain: {}'.format(pixel_gain[y,0]))
            if (write_irs=='y'):
                if correct_pixel=='y':
                    calib_gain = pixel_gain
                elif correct_atlas=='y':
                    calib_gain = atlas_gain
                else:
                    calib_gain = 1.0
                sc.correct_spec(calib_gain, filOUT=f2nam+'_'+str(j))
        
    if (write_irs=='y' and Nmc>1):
        for iobs in trange(Nobs, leave=False,
                           desc='Calculating uncertainties for IRS slits'):
            if mod=='SL':
                f2nam = path_tmp+'SL/'+parobs[iobs][0]
            elif mod=='LL':
                f2nam = path_tmp+'LL/'+parobs[iobs][0]
            else:
                f2nam = None
            fits_Nmc = pathlib.Path(f2nam+'_'+str(Nmc)+fitsext)
            if fits_Nmc.exists():
                mcimage = []
                for j in range(Nmc):
                    ds = read_fits(f2nam+'_'+str(j+1))
                    mcimage.append(ds.data)
                mcimage = np.array(mcimage)
                unc = np.nanstd(mcimage, axis=0)
                write_fits(f2nam+'_unc', ds.header, unc, ds.wave)


##----------------------------------------------------------

##                   Stitch IRC-IRS spectra

##----------------------------------------------------------
# concat_irs = input("Stitch SL-LL spectra of IRS sub-maps (y/n)? ")
concat_irs = None
concat_irc = input("Stitch SL-LL spectra that are reprojected to IRC slits (y/n)? ")
concat_mir = input("Stitch IRC-IRS spectra (y/n)? ")

if concat_irs=='y':
    for j in trange(Nmc+1,# leave=False
                    desc='Stitching SL-LL (IRS frame) [MC]'):
        if j==0:
            concatenate((path_out+src+'_SL',path_out+src+'_LL'),
                        path_out+src+'_IRS', wsort=False,# wrange=wrange,
                        keepfrag=False, cropedge=False)
        else:
            concatenate((path_tmp+src+'_SL_'+str(j),path_tmp+src+'_LL_'+str(j)),
                        path_tmp+src+'_IRS'+str(j), wsort=False,# wrange=wrange,
                        keepfrag=False, cropedge=False)
    
    ## Uncertainty cube
    fits_Nmc = pathlib.Path(path_tmp+src+'_IRS_'+str(Nmc)+fitsext)
    if Nmc>1 and fits_Nmc.exists():
        mcimage = []
        for j in range(Nmc):
            ds = read_fits(path_tmp+src+'_IRS_'+str(j+1))
            mcimage.append(ds.data)
        mcimage = np.array(mcimage)
        unc = np.nanstd(mcimage, axis=0)
        write_fits(path_out+src+'_IRS_unc', ds.header, unc, ds.wave)

    ## Extra stitching of SL1-LL2 (due to MIPS1 map deficit/saturation)
    ds = read_fits(path_out+src+'_IRS', path_out+src+'_IRS_unc')
    Nw, Ny, Nx = ds.data.shape
    for x in range(Nx):
        for y in range(Ny):
            i1 = closest(ds.wave, 14.29, 'left') + 1
            gain = ds.data[i1-1,y,x]/ds.data[i1,y,x]
            ds.data[i1:,y,x] = ds.data[i1:,y,x] * gain
            ds.unc[i1:,y,x] = ds.unc[i1:,y,x] * gain
    write_fits(path_out+src+'_IRS', ds.header, ds.data, ds.wave)
    write_fits(path_out+src+'_IRS_unc', ds.header, ds.unc, ds.wave)

if concat_irc=='y':
    for iobs in trange(Nobs,# leave=False,
                       desc='Stitching SL-LL (IRC slits)'):
        xscale, yscale = read_hdf5(filog+parobs[iobs][0], 'Super pixel size')

        for j in trange(Nmc+1, leave=False,
                        desc=' - '+parobs[iobs][0]+' [MC]'):
            if j==0:
                concatenate((path_tmp+'SL/'+parobs[iobs][0]+'_'+str(j),
                             path_tmp+'LL/'+parobs[iobs][0]+'_'+str(j)),
                            out_irs[iobs], wsort=False,# wrange=wrange,
                            keepfrag=False, cropedge=False)
                igroupixel(out_irs[iobs],
                           xscale=xscale, yscale=yscale,
                           filOUT=out_irs[iobs])
            else:
                concatenate((path_tmp+'SL/'+parobs[iobs][0]+'_'+str(j),
                             path_tmp+'LL/'+parobs[iobs][0]+'_'+str(j)),
                            path_tmp+src+'_'+parobs[iobs][0]+'_IRS_'+str(j),
                            wsort=False,# wrange=wrange,
                            keepfrag=False, cropedge=False)
                igroupixel(path_tmp+src+'_'+parobs[iobs][0]+'_IRS_'+str(j),
                           xscale=xscale, yscale=yscale,
                           filOUT=path_tmp+src+'_'+parobs[iobs][0]+'_IRS_'+str(j))
        
        ## Uncertainty cube
        fits_Nmc = pathlib.Path(path_tmp+src+'_'+parobs[iobs][0]+'_IRS_'+str(Nmc)+fitsext)
        if Nmc>1 and fits_Nmc.exists():
            mcimage = []
            for j in range(Nmc):
                ds = read_fits(path_tmp+src+'_'+parobs[iobs][0]+'_IRS_'+str(j+1))
                mcimage.append(ds.data)
            mcimage = np.array(mcimage)
            unc = np.nanstd(mcimage, axis=0)
            write_fits(out_irs[iobs]+'_unc', ds.header, unc, ds.wave)

        ## Extra stitching of SL1-LL2 (due to MIPS1 map deficit/saturation)
        ds = read_fits(out_irs[iobs], out_irs[iobs]+'_unc')
        Nw, Ny, Nx = ds.data.shape
        for x in range(Nx):
            for y in range(Ny):
                i1 = closest(ds.wave, 14.29, 'left') + 1
                gain = ds.data[i1-1,y,x]/ds.data[i1,y,x]
                ds.data[i1:,y,x] = ds.data[i1:,y,x] * gain
                ds.unc[i1:,y,x] = ds.unc[i1:,y,x] * gain
        write_fits(out_irs[iobs], ds.header, ds.data, ds.wave)
        write_fits(out_irs[iobs]+'_unc', ds.header, ds.unc, ds.wave)


if concat_mir=='y':
    for iobs in trange(Nobs,# leave=False,
                       desc='Stitching IRC-IRS'):
        for j in trange(Nmc+1, leave=False,
                        desc=' - '+parobs[iobs][0]+' [MC]'):
            if j==0:
                concatenate((out_irc[iobs],out_irs[iobs]),
                            path_out+src+'_'+parobs[iobs][0],
                            wsort=False,# wrange=wrange,
                            keepfrag=True, cropedge=False) # keepfrag=True to show IRC only spectra
                ## Added a NaN at 5 micron
                ds = read_fits(path_out+src+'_'+parobs[iobs][0])
                i1 = closest(ds.wave, 5., 'left') + 1
                # i2 = closest(ds.wave, 20.5, 'left') + 1
                i2 = closest(ds.wave, 40, 'left') + 1
                data = np.insert(ds.data, i1, np.nan, axis=0)
                wave = np.insert(ds.wave, i1, 5.)
                write_fits(path_out+src+'_'+parobs[iobs][0],
                           ds.header, data[:i2,:,:], wave[:i2])
            else:
                concatenate((fits_irc[iobs]+'_'+str(j),
                             path_tmp+src+'_'+parobs[iobs][0]+'_IRS_'+str(j)),
                            path_tmp+src+'_'+parobs[iobs][0]+'_'+str(j),
                            wsort=False,# wrange=wrange,
                            keepfrag=True, cropedge=False)
                ## Added a NaN at 5 micron
                ds = read_fits(path_tmp+src+'_'+parobs[iobs][0]+'_'+str(j))
                data = np.insert(ds.data, i1, np.nan, axis=0)
                wave = np.insert(ds.wave, i1, 5.)
                write_fits(path_tmp+src+'_'+parobs[iobs][0]+'_'+str(j),
                           ds.header, data[:i2,:,:], wave[:i2])
        ## Uncertainty cube
        fits_Nmc = pathlib.Path(path_tmp+src+'_'+parobs[iobs][0]+'_'+str(Nmc)+fitsext)
        if Nmc>1 and fits_Nmc.exists():
            mcimage = []
            for j in range(Nmc):
                ds = read_fits(path_tmp+src+'_'+parobs[iobs][0]+'_'+str(j+1))
                mcimage.append(ds.data)
            mcimage = np.array(mcimage)
            unc = np.nanstd(mcimage, axis=0)
            write_fits(path_out+src+'_'+parobs[iobs][0]+'_unc', ds.header, unc, ds.wave)


##----------------------------------------------------------

##                      Plot spectra

##----------------------------------------------------------
plot_mir = input("Plot MIR spectra (y/n)? ")
if plot_mir=='y':

    ## Individual spectra (IRS)
    ##--------------------------
    for iobs in range(Nobs):
        ds = read_fits(path_out+src+'_'+parobs[iobs][0], path_out+src+'_'+parobs[iobs][0]+'_unc')
        Nw, Ny, Nx = ds.data.shape
        iw5 = closest(ds.wave, 5, 'left') + 1
        iw14 = closest(ds.wave, 14.29, 'left') + 1
        iw20 = closest(ds.wave, 20.67, 'left') + 1
        # lolims = np.zeros((Nw,Ny,Nx))
        clib = ['c','m','y','r','orange','b','g']
        xtic = [2, 3, 4, 5, 6, 7, 8, 10, 12, 15, 20]
        xtic_min = np.arange(2.5, 20, .1)
        # xtic = [2, 3, 4, 5, 6, 7, 8, 10, 12, 15, 20, 30, 40]
        # xtic_min = [np.arange(2.5, 20, .1), np.arange(20, 40, 1)]
        # xtic_min = np.concatenate(xtic_min)

        for x in range(Nx):
            for y in range(Ny):
                ## Add calibration error
                ds.unc[:,y,x] = np.sqrt( ds.unc[:,y,x]**2 + (ds.data[:,y,x]*.05)**2 ) # IRS 5% (Carey2010)

                for k in range (Nw):
                    ## Exclude negtive flux
                    if ds.data[k,y,x]<0:
        #                 if iobs==10:
                        ds.data[k,y,x] = -ds.data[k,y,x]
        #                 else:
        #                     ds.data[k,y,x] = np.nan
        #                     ds.unc[k,y,x] = np.nan
        #             ## Set lower limits of unc that are larger than data
        #             if ds.data[k,y,x]<=ds.unc[k,y,x]:
        #                 lolims[k,y,x] = True

                ## Lower limit
                # ds.data[:,y,x][ds.data[:,y,x]<1e-2] = np.nan
                # ds.unc[:,y,x][ds.data[:,y,x]<1e-2] = 0

        xscale, yscale = read_hdf5(filog+parobs[iobs][0], 'Super pixel size')
        for y in range(Ny):
            ys = pix2sup(y, yscale, origin=0)
            subname = parobs[iobs][0][0]+str(ys+1)
            
            ## Remove artifacts (mainly in LL1)
            maskspec = np.full(Nw, True)
            if (subname[0]=='A' and int(subname[1])>2) \
               or (subname[0]=='B' and int(subname[1])>0) \
               or subname=='C2' or subname=='C3' or subname=='C4' \
               or subname=='D3' or subname=='D4' \
               or subname=='F1' \
               or subname[0]=='G' \
               or subname[0]=='H' \
               or subname=='I1':
                maskspec[iw20:] = False
            elif subname[0]=='E' or subname[0]=='N':
                maskspec[:] = False

            if (maskspec.any() and y%yscale==0):
                p0 = pplot(ds.wave[iw5:], ds.data[iw5:,y,0], yerr=ds.unc[iw5:,y,0],
                           xlog=1, ylog=1, nonposx='clip', nonposy='clip',
                           c='k', lw=2, ec='r', label=subname,
                           xlim=(5,21.),
                           xlabel=r'$\rm Wavelengths\ \lambda\ (\mu m)$',
                           ylabel=r'$\rm F_{\nu}\ (MJy/sr)$',
                           # xtk=xtic, xtkmi=xtic_min,
                           xtkform='mylog', ytkform='log_sci',
                           figsize=(12,9), loc='upper left', legendalpha=0,
                           # title=src+'_'+subname,
                           titlesize=20, xysize=20, tksize=20, legendsize=20)

                # p0.ax.text(21, ds.data[iw5,y,0]*1.2, 'LL1 fringe pattern',
                #            size=20, c='k',transform=p0.ax.transData)
                # p0.ax.annotate(xy=(20.67,ds.data[iw5,y,0]*.9), transform=p0.ax.transData,
                #                xytext=(38,ds.data[iw5,y,0]*.9),
                #                text='', c='k', arrowprops=dict(arrowstyle='<->'))
                
                p0.save(path_fig+'IRS_'+subname, transparent=True, figtight=True)

    ## Overlapped in one figure
    ##--------------------------
    for iobs in range(Nobs):
        ds = read_fits(path_out+src+'_'+parobs[iobs][0], path_out+src+'_'+parobs[iobs][0]+'_unc')
        Nw, Ny, Nx = ds.data.shape
        iw5 = closest(ds.wave, 5, 'left') + 1
        iw14 = closest(ds.wave, 14.29, 'left') + 1
        iw20 = closest(ds.wave, 20.67, 'left') + 1
        # lolims = np.zeros((Nw,Ny,Nx))
        clib = ['c','m','y','r','orange','b','g']
        xtic = [2, 3, 4, 5, 6, 7, 8, 10, 12, 15, 20]
        xtic_min = np.arange(2.5, 20, .1)
        # xtic = [2, 3, 4, 5, 6, 7, 8, 10, 12, 15, 20, 30, 40]
        # xtic_min = [np.arange(2.5, 20, .1), np.arange(20, 40, 1)]
        # xtic_min = np.concatenate(xtic_min)

        for x in range(Nx):
            for y in range(Ny):
                ## Add calibration error
                ds.unc[:,y,x] = np.sqrt( ds.unc[:,y,x]**2 + (ds.data[:,y,x]*.05)**2 ) # IRS 5% (Carey2010)

                for k in range (Nw):
                    ## Exclude negtive flux
                    if ds.data[k,y,x]<0:
        #                 if iobs==10:
                        ds.data[k,y,x] = -ds.data[k,y,x]
        #                 else:
        #                     ds.data[k,y,x] = np.nan
        #                     ds.unc[k,y,x] = np.nan
        #             ## Set lower limits of unc that are larger than data
        #             if ds.data[k,y,x]<=ds.unc[k,y,x]:
        #                 lolims[k,y,x] = True

                ## Lower limit
                # ds.data[:,y,x][ds.data[:,y,x]<1e-2] = np.nan
                # ds.unc[:,y,x][ds.data[:,y,x]<1e-2] = 0

        xscale, yscale = read_hdf5(filog+parobs[iobs][0], 'Super pixel size')
        for y in range(Ny):
            ys = pix2sup(y, yscale, origin=0)
            subname = parobs[iobs][0][0]+str(ys+1)
            
            ## Remove artifacts (mainly in LL1)
            maskspec = np.full(Nw, True)
            if (subname[0]=='A' and int(subname[1])>2) \
               or (subname[0]=='B' and int(subname[1])>0) \
               or subname=='C2' or subname=='C3' or subname=='C4' \
               or subname=='D3' or subname=='D4' \
               or subname=='F1' \
               or subname[0]=='G' \
               or subname[0]=='H' \
               or subname=='I1':
                maskspec[iw20:] = False
            elif subname[0]=='E' or subname[0]=='N':
                maskspec[:] = False

            if (maskspec.any() and y%yscale==0):
                xlim = (2.5,21.)
                # xlim = (2,40)
                if subname[0]=='A':
                    ylim = (1e-1,5e6)
                elif subname[0]=='B':
                    ylim = (1e-1,3e5)
                elif subname[0]=='C':
                    ylim = (1e-1,1e6)
                elif subname[0]=='D':
                    ylim = (1e-2,1e4)
                else:
                    ylim = (None,None)
                
                if ys==0:
                    p = pplot(xlog=1, ylog=1, nonposx='clip', nonposy='clip',
                              lw=2, ec='grey', elw=1,
                              xlim=xlim, ylim=ylim,# capsize=2,
                              xlabel=r'$\rm Wavelengths\ \lambda\ (\mu m)$',
                              ylabel=r'$\rm F_{\nu}\ (MJy/sr)$',
                              xtk=xtic, xtkmi=xtic_min, xtkform='mylog', ytkform='log_sci',
                              figsize=(15,8), loc='upper left', legendalpha=0,
                              # title=src+' '+parobs[iobs][0]+' spectra',
                              titlesize=20, xysize=20, tksize=20, legendsize=20)

                    # if iobs==0:
                    #     p.ax.text(.85,.05,'(a)',size=30,c='grey',
                    #               transform=p.ax.transAxes) # for the use of Hu_thesis
                    # elif iobs==6:
                    #     p.ax.text(.85,.05,'(b)',size=30,c='grey',
                    #               transform=p.ax.transAxes) # for the use of Hu_thesis
                    
                if subname[1]=='1' or subname[1]=='7' \
                   or subname=='B5' or subname=='C5' or subname=='C6' \
                   or subname=='D3':
                    ds.data[:,y,0] *= .1
                    ds.unc[:,y,0] *= .1
                    subname += '*0.1'
                elif subname=='A6':
                    ds.data[:,y,0] *= .01
                    ds.unc[:,y,0] *= .01
                    subname += '*0.01'
                elif subname=='B6':
                    ds.data[:,y,0] *= .08
                    ds.unc[:,y,0] *= .08
                    subname += '*0.08'
                # else:
                #     ds.data[:,y,0] *= 1
                #     ds.unc[:,y,0] *= 1
                #     subname += ''

                p.add_plot(ds.wave[maskspec], ds.data[:,y,0][maskspec],
                           yerr=ds.unc[:,y,0][maskspec],# lolims=lolims[:,y,0],
                           c=clib[ys], lw=2, ec='grey', elw=1, label=subname)

        if (subname[0]!='E' and subname[0]!='N'):
            p.ax.text(p.transData2Axes((3,1))[0], 0.07, 'IRC',
                      size=20, c='k', transform=p.ax.transAxes)
            p.ax.text(p.transData2Axes((6,1))[0], 0.07, 'SL2',
                      size=20, c='k', transform=p.ax.transAxes)
            p.ax.text(p.transData2Axes((9,1))[0], 0.07, 'SL1',
                      size=20, c='k', transform=p.ax.transAxes)
            p.ax.text(p.transData2Axes((16,1))[0], 0.07, 'LL2',
                      size=20, c='k', transform=p.ax.transAxes)
            # p.ax.text(p.transData2Axes((25,1))[0], 0.07, 'LL1',
            #           size=20, c='k', transform=p.ax.transAxes)
            p.ax.annotate(xy=(p.transData2Axes((2.5,1))[0],0.05),
                          xycoords=p.ax.transAxes,
                          xytext=(p.transData2Axes((5,1))[0],0.05),
                          text='', c='k', arrowprops=dict(arrowstyle='<->'))
            p.ax.annotate(xy=(p.transData2Axes((5.21,1))[0],0.05),
                          xycoords=p.ax.transAxes,
                          xytext=(p.transData2Axes((7.56,1))[0],0.05),
                          text='', c='k', arrowprops=dict(arrowstyle='<->'))
            p.ax.annotate(xy=(p.transData2Axes((7.57,1))[0],0.05),
                          xycoords=p.ax.transAxes,
                          xytext=(p.transData2Axes((14.28,1))[0],0.05),
                          text='', c='k', arrowprops=dict(arrowstyle='<->'))
            p.ax.annotate(xy=(p.transData2Axes((14.29,1))[0],0.05),
                          xycoords=p.ax.transAxes,
                          xytext=(p.transData2Axes((20.66,1))[0],0.05),
                          text='', c='k', arrowprops=dict(arrowstyle='<->'))
            # p.ax.annotate(xy=(p.transData2Axes((20.67,1))[0],0.05),
            #               xycoords=p.ax.transAxes,
            #               xytext=(p.transData2Axes((38.0,1))[0],0.05),
            #               text='', c='k', arrowprops=dict(arrowstyle='<->'))
        
            p.save(path_fig+src+'_'+subname[0], transparent=True, figtight=True)

