#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""

Build stitched Spitzer/IRS spectral cube
Inter-calibrate Spitzer/IRS spectra with Spitzer/IRAC4 (SL) & MIPS1 (LL)

"""

# import logging, sys
# logging.basicConfig(stream=sys.stderr, level=logging.ERROR)
# print(logging.getLogger())
# logging.disable(sys.maxsize) # disable IDL print
import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning) 

from tqdm import tqdm, trange

import os
import numpy as np
from scipy.optimize import curve_fit
from scipy.integrate import trapz
from matplotlib.ticker import ScalarFormatter, NullFormatter

## laputan
from laputan.inout import fclean, fitsext, read_fits, write_fits
from laputan.imaging import ( iconvolve, iswarp, iuncert, concatenate,
                              imontage, improve, Jy_per_pix_to_MJy_per_sr )
from laputan.astrom import fixwcs
from laputan.arrays import closest
from laputan.calib import intercalib
from laputan.maths import f_lin, f_lin0
from laputan.plots import pplot


##---------------------------
##       Preliminaries
##---------------------------
## Local
from buildinfo import ( src, Nmc, verbose, coadd_tool,
                        chnl, fits_irs, out_irs, out_irc, 
                        path_idl, path_ker, path_conv, path_phot, path_cal,
                        fits_ker, csv_ker, path_tmp, path_fig
)
Nch = len(chnl)
Nmc = 20

phot = ['IRAC4', 'MIPS1']
Nphot = len(phot)

##---------------------------
##    Coadd observations
##---------------------------
swp_coadd = input("Coadd IRS observations (y/n)? ")
if swp_coadd=='y':
    refheader = fixwcs(out_irc+fitsext).header
    swp = iswarp(refheader=refheader,
    # swp = iswarp(sum(fits_irs, []), refheader=refheader,
                 tmpdir=path_tmp, verbose=verbose)
    
    ## Add MC unc
    ##------------
    for i in trange(Nch, #leave=False,
                    desc='<iswarp> IRS Coadding ({} chnl)'.format(Nch)):
        for j in trange(Nmc+1, leave=False,
                        desc='<iswarp> IRS Coadding [MC]'):
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
    # for i in trange(Nch,
    #                 desc='IRS Cal unc ({} chnl)'.format(Nch)):
    #     mcimage = []
    #     for j in trange(Nmc+1, leave=False,
    #                     desc='IRS Reading [MC]'):
    #         if j==0:
    #             hd0 = read_fits(path_tmp+src+'_'+chnl[i]+'_conv')
    #             header = hd0.header
    #             wvl = hd0.wave
    #         else:
    #             hd = read_fits(path_tmp+src+'_'+chnl[i]+'_'+str(j)+'_conv')
    #             mcimage.append(hd.data)
    #     if Nmc>1:
    #         mcimage = np.array(mcimage)
    #         unc = np.nanstd(mcimage, axis=0)
    #         write_fits(path_tmp+src+'_'+chnl[i]+'_conv_unc', header, unc, wvl)

##---------------------------
## Match SL2-SL1 and LL2-LL1
##---------------------------
do_match = input("Match SL2-SL1 and LL2-LL1 (y/n)? ")
if do_match=='y':
    ## SL3
    data_sl3 = read_fits(path_tmp+src+'_SL3_conv').data
    data_sl2 = read_fits(path_tmp+src+'_SL2_conv').data
    data_sl1 = read_fits(path_tmp+src+'_SL1_conv').data
    wvl_sl3 = read_fits(path_tmp+src+'_SL3_conv').wave
    wvl_sl2 = read_fits(path_tmp+src+'_SL2_conv').wave
    wvl_sl1 = read_fits(path_tmp+src+'_SL1_conv').wave
    iwmax_sl2 = closest(wvl_sl3, wvl_sl2[-1])
    iwmin_sl1 = closest(wvl_sl3, wvl_sl1[0])
    left_sl3 = trapz(data_sl3[:iwmax_sl2], wvl_sl3[:iwmax_sl2],
                     dx=wvl_sl3[1]-wvl_sl3[0], axis=0)
    right_sl3 = trapz(data_sl3[iwmin_sl1:], wvl_sl3[iwmin_sl1:],
                      dx=wvl_sl3[1]-wvl_sl3[0], axis=0)
    ## LL3
    data_ll3 = read_fits(path_tmp+src+'_LL3_conv').data
    data_ll2 = read_fits(path_tmp+src+'_LL2_conv').data
    data_ll1 = read_fits(path_tmp+src+'_LL1_conv').data
    wvl_ll3 = read_fits(path_tmp+src+'_LL3_conv').wave
    wvl_ll2 = read_fits(path_tmp+src+'_LL2_conv').wave
    wvl_ll1 = read_fits(path_tmp+src+'_LL1_conv').wave
    iwmax_ll2 = closest(wvl_ll3, wvl_ll2[-1])
    iwmin_ll1 = closest(wvl_ll3, wvl_ll1[0])
    left_ll3 = trapz(data_ll3[:iwmax_ll2], wvl_ll3[:iwmax_ll2],
                     dx=wvl_ll3[1]-wvl_ll3[0], axis=0)
    right_ll3 = trapz(data_ll3[iwmin_ll1:], wvl_ll3[iwmin_ll1:],
                      dx=wvl_ll3[1]-wvl_ll3[0], axis=0)
    
    iwmin_sl3 = closest(wvl_sl2, wvl_sl3[0]) # SL2 index
    iwmax_sl3 = closest(wvl_sl1, wvl_sl3[-1]) # SL1 index
    iwmin_ll3 = closest(wvl_ll2, wvl_ll3[0]) # LL2 index
    iwmax_ll3 = closest(wvl_ll1, wvl_ll3[-1]) # LL1 index
    right_sl2 = trapz(data_sl2[iwmin_sl3:], wvl_sl2[iwmin_sl3:],
                      dx=wvl_sl2[1]-wvl_sl2[0], axis=0)
    left_sl1 = trapz(data_sl1[:iwmax_sl3], wvl_sl1[:iwmax_sl3],
                     dx=wvl_sl1[1]-wvl_sl1[0], axis=0)
    right_ll2 = trapz(data_ll2[iwmin_ll3:], wvl_ll2[iwmin_ll3:],
                      dx=wvl_ll2[1]-wvl_ll2[0], axis=0)
    left_ll1 = trapz(data_ll1[:iwmax_ll3], wvl_ll1[:iwmax_ll3],
                     dx=wvl_ll1[1]-wvl_ll1[0], axis=0)

    Ny, Nx = data_sl2.shape[1:]
    offset = np.zeros((4,Ny,Nx))
    offset[0] = left_sl3 - right_sl2
    offset[1] = right_sl3 - left_sl1
    offset[2] = left_ll3 - right_ll2
    offset[3] = right_ll3 - left_ll1
            
    ## Match SL2 to SL3
    intercalib(path_tmp+src+'_SL2_conv').specorrect(filOUT=path_tmp+src+'_SL2_off',
                                                    offset=offset[0])
    ## Match SL1 to SL3
    intercalib(path_tmp+src+'_SL1_conv').specorrect(filOUT=path_tmp+src+'_SL1_off',
                                                    offset=offset[1])
    ## Match LL2 to LL3
    intercalib(path_tmp+src+'_LL2_conv').specorrect(filOUT=path_tmp+src+'_LL2_off',
                                                    offset=offset[2])
    ## Match LL1 to LL3
    intercalib(path_tmp+src+'_LL1_conv').specorrect(filOUT=path_tmp+src+'_LL1_off',
                                                    offset=offset[3])
##---------------------------
##     Concatenate IRS
##---------------------------
concat_irs = input("Stitch IRS spectra (y/n)? ")
if concat_irs=='y':

    # ## SL3
    # data_sl3 = read_fits(path_tmp+src+'_SL3_conv').data
    # wvl_sl3 = read_fits(path_tmp+src+'_SL3_conv').wave
    # ind_sl = closest(wvl_sl3, 7.56)
    # left_sl3 = trapz(data_sl3[:ind_sl], wvl_sl3[:ind_sl],
    #                  dx=wvl_sl3[1]-wvl_sl3[0], axis=0)
    # right_sl3 = trapz(data_sl3[ind_sl:], wvl_sl3[ind_sl:],
    #                   dx=wvl_sl3[1]-wvl_sl3[0], axis=0)
    # ## LL3
    # data_ll3 = read_fits(path_tmp+src+'_LL3_conv').data
    # wvl_ll3 = read_fits(path_tmp+src+'_LL3_conv').wave
    # ind_ll = closest(wvl_ll3, 20.66)
    # left_ll3 = trapz(data_ll3[:ind_ll], wvl_ll3[:ind_ll],
    #                  dx=wvl_ll3[1]-wvl_ll3[0], axis=0)
    # right_ll3 = trapz(data_ll3[ind_ll:], wvl_ll3[ind_ll:],
    #                   dx=wvl_ll3[1]-wvl_ll3[0], axis=0)
    
    for j in trange(Nmc+1,
                    desc='IRS concatenation [MC]'):
        files = []
        if j==0:
            for i in range(Nch-2):
                files.append(path_tmp+src+'_'+chnl[i]+'_conv')
    
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
            # Nw, Ny, Nx = data.shape
            header = hd.header
            wvl = hd.wave
            # offset = np.zeros((4,Ny,Nx))
        else:
            for i in range(Nch-2):
                files.append(path_tmp+src+'_'+chnl[i]+'_'+str(j)+'_conv')
    
            concatenate(files, out_irs+'_'+str(j), wsort=False)
            data = read_fits(out_irs+'_'+str(j)).data

        # iwi_sl3 = closest(wvl, wvl_sl3[0])
        # iw_sl = closest(wvl, 7.56)
        # iws_sl3 = closest(wvl, wvl_sl3[-1])
        # iwi_ll3 = closest(wvl, wvl_ll3[0])
        # iw_ll = closest(wvl, 20.66)
        # iws_ll3 = closest(wvl, wvl_ll3[-1])
        # dw = wvl[1]-wvl[0]
        # right_sl2 = trapz(data[iwi_sl3:iw_sl], wvl[iwi_sl3:iw_sl], dx=dw, axis=0)
        # left_sl1 = trapz(data[iw_sl:iws_sl3], wvl[iw_sl:iws_sl3], dx=dw, axis=0)
        # right_ll2 = trapz(data[iwi_ll3:iw_ll], wvl[iwi_ll3:iw_ll], dx=dw, axis=0)
        # left_ll1 = trapz(data[iw_ll:iws_ll3], wvl[iw_ll:iws_ll3], dx=dw, axis=0)
        # offset[0] = left_sl3 - right_sl2
        # offset[1] = right_sl3 - left_sl1
        # offset[2] = left_ll3 - right_ll2
        # offset[3] = right_ll3 - left_ll1
                
        # write_fits(out_irs+'_'+str(j), header, data, wvl)
        # ## Match SL2 to SL3
        # intercalib(out_irs+'_'+str(j)).specorrect(filOUT=out_irs+'_'+str(j),
        #                                           offset=offset[0],
        #                                           wlim=(5.21, 7.56))
        # ## Match SL1 to SL3
        # intercalib(out_irs+'_'+str(j)).specorrect(filOUT=out_irs+'_'+str(j),
        #                                           offset=offset[1],
        #                                           wlim=(7.57, 14.28))
        # ## Match LL2 to LL3
        # intercalib(out_irs+'_'+str(j)).specorrect(filOUT=out_irs+'_'+str(j),
        #                                           offset=offset[2],
        #                                           wlim=(14.29, 20.66))
        # ## Match LL1 to LL3
        # intercalib(out_irs+'_'+str(j)).specorrect(filOUT=out_irs+'_'+str(j),
        #                                           offset=offset[3],
        #                                           wlim=(20.67, 38.00))
    
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

##---------------------------
##    Prepare photometry
##---------------------------
pre_calib1 = input("Run inter-calibration prepipeline (DustPedia) (y/n)? ")
pre_calib2 = input("Run inter-calibration prepipeline (SINGS) (y/n)? ")
int_calib = input("Run inter-calibration (y/n)? ")

if not os.path.exists(path_tmp+'calib/'):
    os.makedirs(path_tmp+'calib/')
    
for i in range(Nphot):
    
    if phot[i]=='IRAC4':
        phot_ker = path_ker+'Kernel_HiRes_IRAC_8.0_to_Gauss_06.0'
    elif phot[i]=='MIPS1':
        phot_ker = path_ker+'Kernel_HiRes_MIPS_24_to_Gauss_06.0'
    
    raw_phot1 = path_phot+src+'_'+phot[i]+'_DP'
    fits_phot1 = path_cal+src+'_'+phot[i]+'_DP'
    tmp_phot1 = path_tmp+'calib/'+src+'_'+phot[i]+'_DP'
    
    raw_phot2 = path_phot+src+'_'+phot[i]+'_SINGS'
    wgt_phot2 = path_phot+src+'_'+phot[i]+'_SINGS_wgt'
    fits_phot2 = path_cal+src+'_'+phot[i]+'_SINGS'
    tmp_phot2 = path_tmp+'calib/'+src+'_'+phot[i]+'_SINGS'
    
    fits_spec = path_cal+src+'_'+phot[i]+'_IRS'
    tmp_spec = path_tmp+'calib/'+src+'_'+phot[i]+'_IRS'
    
    ## DustPedia (phot1)
    ##===================
    if pre_calib1=='y':

        if i==0:
            ## Convert phot unit (Jy/pix -> MJy/sr)
            ## This step should be before iswarp reprojection (pixscale changed)
            Jy_per_pix_to_MJy_per_sr(raw_phot1, filOUT=path_tmp+src+'_'+phot[i]+'_DP')
            Jy_per_pix_to_MJy_per_sr(raw_phot1+'_unc',
                                     filOUT=path_tmp+src+'_'+phot[i]+'_DP_unc')
        
            ## Reproject phot
            refheader = fixwcs(out_irc+fitsext).header
            mtg = imontage('exact', tmpdir=path_tmp)
            mtg.reproject_mc(path_tmp+src+'_'+phot[i]+'_DP', refheader=refheader,
                             dist='norm', Nmc=Nmc, filOUT=tmp_phot1)
            
            ## Convolve phot
            for j in trange(Nmc+1, leave=False,
                            desc=phot[i]+': DP conv [MC]'):
                if j==0:
                    conv = iconvolve(tmp_phot1, phot_ker, csv_ker, filOUT=fits_phot1)
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
    if pre_calib2=='y':
        
        ## Create uncertainty via weight map (suppose uniform contribution)
        ## The weight maps contain the information on the number of frames 
        ## that were used to create the science mosaics at each pixel (value= # frames x10); 
        ## the pixel size of the weight maps is the same as the science mosaics.
        ## -- SINGS v5 release doc
        # if phot[i]== 'IRAC4':
        #     bg = read_fits(raw_phot2).data[502:542,1067:1107]
        #     bg_wgt = read_fits(wgt_phot2).data[502:542,1067:1107] / 10.
        # elif phot[i]=='MIPS1':
        #     bg = read_fits(raw_phot2).data[450:490,1606:1646]
        #     bg_wgt = read_fits(wgt_phot2).data[450:490,1606:1646] / 10.
        # iuncert(raw_phot2, filOUT=raw_phot2+'_unc',
        #         filWGT=wgt_phot2, wfac=.1,
        #         BG_image=bg, BG_weight=bg_wgt)

        ## Reproject phot
        refheader = fixwcs(out_irc+fitsext).header
        mtg = imontage('exact', tmpdir=path_tmp)
        mtg.reproject_mc(raw_phot2, refheader=refheader,
                         dist='norm', Nmc=Nmc,
                         filOUT=tmp_phot2)
    
        ## Convolve phot
        for j in trange(Nmc+1,# leave=False,
                        desc=phot[i]+': SINGS conv [MC]'):
            if j==0:
                conv = iconvolve(tmp_phot2, phot_ker, csv_ker, filOUT=fits_phot2)
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
            # print('Calculating uncertainty cube...')
            mcimage = np.array(mcimage)
            unc = np.nanstd(mcimage, axis=0)
            write_fits(fits_spec+'_unc', ic.hdr, unc)

##---------------------------
##     Inter-Calibration
##---------------------------

    ## Data
    ##------
    if i==0:
        data_phot1 = read_fits(fits_phot1).data
        unc_phot1 = read_fits(fits_phot1+'_unc').data
        pix_phot1 = data_phot1[:,:].reshape((-1,))
        pix_phot1_unc = unc_phot1[:,:].reshape((-1,))
        
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
            
    if int_calib=='y':
        
        ## Plot
        ##------
        pix_phot2 = data_phot2[:,:].reshape((-1,))
        pix_spec = data_spec[:,:].reshape((-1,))
        pix_phot2_unc = unc_phot2[:,:].reshape((-1,))
        pix_spec_unc = unc_spec[:,:].reshape((-1,))
        
        # xgrid = np.arange(0,1e3,1)
        xgrid = np.logspace(-6,4,10000)
        
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
            p0 = pplot(pix_phot1,  pix_phot2,
                       yerr=pix_phot2_unc, xerr=pix_phot1_unc,
                       fmt='.', c='k', ec='r', elw=1,
                       xlog=1, ylog=1,
                       xlim=(0,1e4), ylim=(0,1e4),
                       xlab='DustPedia (MJy/sr)', ylab='SINGS (MJy/sr)',
                       figsize=(8,8), title=src+'_'+phot[i])
            p0.save(path_cal+'SINGS-DP_'+phot[i])
            
        ## IRS - SINGS
        p = pplot(pix_phot2, pix_spec,
                  yerr=pix_spec_unc, xerr=pix_phot2_unc,
                  fmt='.', c='k', ec='r', legend='upper left',
                  xlog=1, ylog=1, nonposx='clip', nonposy='clip',
                  xlim=(1e-6,1e4), ylim=(1e-6,1e4),
                  xlab='SINGS (MJy/sr)', ylab='IRS (MJy/sr)',
                  figsize=(8,8), title=src+'_'+phot[i])
        
        ## Linear fit
        if i==0:
            popt, pcov = curve_fit(f_lin, pix_spec[mask], pix_phot2[mask])
                                   # sigma=pix_phot2_unc[mask])
            calib_fac = popt[0]
            calib_off = popt[1]
            
            # popt, pcov = curve_fit(f_lin, pix_phot2[mask], pix_spec[mask],
            #                        sigma=pix_spec_unc[mask])
            # calib_fac = 1. / popt[0]
            # calib_off = -popt[1] / popt[0]
            
            print('Inter-Calibration ('+phot[i]+') factor = {:.4}'.format(calib_fac))
            print('Inter-Calibration ('+phot[i]+') offset = {:.4}'.format(calib_off))
            # label = phot[i]+' calib fac = {:.4}'.format(calib_fac)
            label = phot[i]+': y={0:.4}x+{1:.4}'.format(calib_fac, calib_off)
            # label = phot[i]+': y={0:.4}x+{1:.4}'.format(*popt)
            
            p.add_plot(f_lin(xgrid, *popt), xgrid,
                       c='y', ls='-', label=label)
        else:
            mask_mips1 = ~np.ma.masked_where(pix_phot2<1.e0, pix_phot2).mask
            mask = np.ma.array(mask,
                               mask=np.logical_and(mask,mask_mips1)).mask
            # popt, pcov = curve_fit(f_lin, pix_spec[mask], pix_phot2[mask],
            #                        sigma=pix_phot2_unc[mask])
            # calib_fac = popt[0]
            # calib_off = popt[1]

            popt, pcov = curve_fit(f_lin, pix_phot2[mask], pix_spec[mask],
                                   sigma=pix_spec_unc[mask])
            calib_fac = 1. / popt[0]
            calib_off = -popt[1] / popt[0]
            
            print('Inter-Calibration ('+phot[i]+') factor = {:.4}'.format(calib_fac))
            print('Inter-Calibration ('+phot[i]+') offset = {:.4}'.format(calib_off))
            # label = phot[i]+' calib fac = {:.4}'.format(calib_fac)
            label = phot[i]+': y={0:.4}x+{1:.4}'.format(calib_fac, calib_off)
            # label = phot[i]+': y={0:.4}x+{1:.4}'.format(*popt)
            
            p.add_plot(xgrid, f_lin(xgrid, *popt),
                       c='y', ls='-', label=label)
        
        p.set_font()
        
        p.save(path_cal+'intercalib_'+phot[i])
        
        
        ## Spectral correction
        ##---------------------
        sc0 = intercalib(out_irs+'_0')
        sc0.specorrect(filOUT=out_irs)
        if i==0:
            wmin = 5.21
            wmax = 14.28
        elif i==1:
            wmin = 14.29
            wmax = 38.00
            
        do_correct = input("Do IRS spectral correction [{}] (y/n)? ".format(phot[i]))
        if do_correct=='y':
            mask_faint = ~np.ma.masked_where((data_phot2-calib_off)/data_spec<1., data_phot2).mask
            
            # calib_fac = data_phot2/data_spec
            imp = improve(out_irs)
            iwmin = closest(imp.wvl, wmin)
            iwmax = closest(imp.wvl, wmax)+1
            im = imp.im
            for iw in range(iwmax-iwmin):
                im[iw+iwmin][mask_faint] = imp.im[iw+iwmin][mask_faint] * calib_fac + calib_off
            write_fits(out_irs, imp.hdr, im, imp.wvl, imp.wmod)

            # sc = intercalib(out_irs)
            # sc.specorrect(filOUT=out_irs,
            #               factor=calib_fac, offset=calib_off,
            #               wlim=(wmin,wmax))

##---------------------------
##       Plot spectra
##---------------------------
do_plot = input("Plot IRS spectra (y/n)? ")
if do_plot=='y':
    ds = read_fits(out_irs, out_irs+'_unc')
    # ds.data = respect().smooth(out_irs, out_irs, lim_unc=1.e2, cmin=2)
    ds0 = read_fits(out_irs+'_0')
    ds_irc = read_fits(out_irc)
    wcen = intercalib().wcenter(phot)
    Nw, Ny, Nx = ds.data.shape
    for x in range(Nx):
        for y in range(Ny):
            if ~np.isnan(ds.data[:,y,x]).any():
                if ~np.isnan(ds_irc.data[0,y,x]).any():
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
    
                    xtic = []
                    xtic_min = np.arange(5., 41., 2.)
                    p.ax.set_xticks(xtic, minor=False) # major
                    p.ax.set_xticks(xtic_min, minor=True) # minor
                    p.ax.xaxis.set_major_formatter(ScalarFormatter()) # major
                    p.ax.xaxis.set_minor_formatter(NullFormatter()) # minor
                    
                    p.save(path_fig+'IRS_('+str(x)+','+str(y)+')')
    
