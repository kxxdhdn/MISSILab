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
from laputan.maths import f_lin, f_lin0, f_lin1
from laputan.plots import pplot


##----------------------------------------------------------

##                     Preliminaries

##----------------------------------------------------------

## Local
from buildinfo import ( src, Nmc, verbose, coadd_tool, coadd_fp,
                        chnl, fits_irs, out_irs, out_irc, 
                        path_idl, path_ker, path_conv, path_phot, path_cal,
                        fits_ker, csv_ker, path_tmp, path_fig
)


## Banner
print('\n============================================================\n')

print('        MIRAGE - Spitzer/IRS cube builder - '+src)

print('\n============================================================\n')

Nch = len(chnl)
## Do not stitch SL3 and LL3
if ('SL3' in chnl) and ('LL3' in chnl):
    Nch_s = Nch - 2
else:
    Nch_s = Nch
    
# Nmc = 2

if ('SH' in chnl) or ('LH' in chnl):
    phot = ['MIPS1']
else:
    phot = ['IRAC4', 'MIPS1']
Nphot = len(phot)

##----------------------------------------------------------

##                   Coadd observations

##----------------------------------------------------------
irs_coadd = input("Coadd IRS observations (y/n)? ")
if irs_coadd=='y':
    
    ## Coadd with adding MC unc
    ##--------------------------
    if coadd_tool=='swarp':
        
        ## <iswarp> coadding
        ##===================
        swp = iswarp(refheader=coadd_fp,
                     tmpdir=path_tmp, verbose=verbose)
        
        for i in trange(Nch, #leave=False,
                        desc='<iswarp> IRS Coadding ({} chnl)'.format(Nch)):
            for j in trange(Nmc+1, leave=False,
                            desc='<iswarp> IRS Coadding [MC]'):
                if j==0:
                    swp.combine(fits_irs[i], combtype='wgt_avg',
                                keepedge=True, cropedge=False,
                                tmpdir=path_tmp+'MC_no/',
                                filOUT=path_tmp+src+'_'+chnl[i])
                else:
                    swp.combine(fits_irs[i], combtype='wgt_avg', 
                                keepedge=True, cropedge=False, uncpdf='norm',
                                tmpdir=path_tmp+'MC_'+str(j)+'/',
                                filOUT=path_tmp+src+'_'+chnl[i]+'_'+str(j))

    elif coadd_tool=='reproject':
                    
        ## <imontage> coadding
        ##=====================
        mtg = imontage('exact', tmpdir=path_tmp, verbose=verbose)
        
        for i in trange(Nch, #leave=False,
                        desc='<imontage> IRS Coadding ({} chnl)'.format(Nch)):
            mtg.coadd(fits_irs[i], refheader=coadd_fp,
                      dist='norm', Nmc=Nmc,
                      filOUT=path_tmp+src+'_'+chnl[i])
    
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
    #             hd = read_fits(path_tmp+src+'_'+chnl[i]+'_conv')
    #             header = hd.header
    #             wvl = hd.wave
    #         else:
    #             hd = read_fits(path_tmp+src+'_'+chnl[i]+'_'+str(j)+'_conv')
    #             mcimage.append(hd.data)
    #     if Nmc>1:
    #         mcimage = np.array(mcimage)
    #         unc = np.nanstd(mcimage, axis=0)
    #         write_fits(path_tmp+src+'_'+chnl[i]+'_conv_unc', header, unc, wvl)

##----------------------------------------------------------

##                   Stitch IRS spectra

##----------------------------------------------------------
concat_irs = input("Stitch IRS spectra (y/n)? ")
lores_match = input(" - Match SL2-SL1 and LL2-LL1 (y/n)? ")
hires_match = input(" - Match SH-LH (y/n)? ")
keep_frag = input(" - Keep fragmentary spectra (y/n)? ")

if keep_frag=='y':
    keepfrag = True
    cropedge = False
else:
    keepfrag = False
    crop_edge_frag = input(" - Crop edge if not keeping frag spectra (y/n)? ")
    if crop_edge_frag=='y':
        cropedge = True
    else:
        cropedge = False

## Match SL2-SL1 and LL2-LL1
##---------------------------
if lores_match=='y':
    ## SL3
    data_sl3 = read_fits(path_tmp+src+'_SL3_conv').data
    data_sl2 = read_fits(path_tmp+src+'_SL2_conv').data
    data_sl1 = read_fits(path_tmp+src+'_SL1_conv').data
    wvl_sl3 = read_fits(path_tmp+src+'_SL3_conv').wave
    wvl_sl2 = read_fits(path_tmp+src+'_SL2_conv').wave
    wvl_sl1 = read_fits(path_tmp+src+'_SL1_conv').wave
    iwmax_sl2 = closest(wvl_sl3, wvl_sl2[-1], side='left') + 1 # SL3 index
    iwmin_sl1 = closest(wvl_sl3, wvl_sl1[0], side='right') # SL3 index
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
    iwmax_ll2 = closest(wvl_ll3, wvl_ll2[-1], side='left') + 1 # LL3 index
    iwmin_ll1 = closest(wvl_ll3, wvl_ll1[0], side='right') # LL3 index
    left_ll3 = trapz(data_ll3[:iwmax_ll2], wvl_ll3[:iwmax_ll2],
                     dx=wvl_ll3[1]-wvl_ll3[0], axis=0)
    right_ll3 = trapz(data_ll3[iwmin_ll1:], wvl_ll3[iwmin_ll1:],
                      dx=wvl_ll3[1]-wvl_ll3[0], axis=0)
    
    iwmin_sl3 = closest(wvl_sl2, wvl_sl3[0], side='right') # SL2 index
    iwmax_sl3 = closest(wvl_sl1, wvl_sl3[-1], side='left') + 1 # SL1 index
    iwmin_ll3 = closest(wvl_ll2, wvl_ll3[0], side='right') # LL2 index
    iwmax_ll3 = closest(wvl_ll1, wvl_ll3[-1], side='left') + 1 # LL1 index
    right_sl2 = trapz(data_sl2[iwmin_sl3:], wvl_sl2[iwmin_sl3:],
                      dx=wvl_sl2[1]-wvl_sl2[0], axis=0)
    left_sl1 = trapz(data_sl1[:iwmax_sl3], wvl_sl1[:iwmax_sl3],
                     dx=wvl_sl1[1]-wvl_sl1[0], axis=0)
    right_ll2 = trapz(data_ll2[iwmin_ll3:], wvl_ll2[iwmin_ll3:],
                      dx=wvl_ll2[1]-wvl_ll2[0], axis=0)
    left_ll1 = trapz(data_ll1[:iwmax_ll3], wvl_ll1[:iwmax_ll3],
                     dx=wvl_ll1[1]-wvl_ll1[0], axis=0)

    Ny, Nx = data_sl2.shape[1:]
    scale = np.zeros((4,Ny,Nx))
    scale[0] = left_sl3 / right_sl2
    scale[1] = right_sl3 / left_sl1
    scale[2] = left_ll3 / right_ll2
    scale[3] = right_ll3 / left_ll1
    ## Display scaling factor map
    mask_sca = ~np.isnan(scale[0])
    print('SL2 to SL3 scaling factor: ', scale[0][mask_sca])
    mask_sca = ~np.isnan(scale[1])
    print('SL1 to SL3 scaling factor: ', scale[1][mask_sca])
    mask_sca = ~np.isnan(scale[2])
    print('LL2 to LL3 scaling factor: ', scale[2][mask_sca])
    mask_sca = ~np.isnan(scale[3])
    print('LL1 to LL3 scaling factor: ', scale[3][mask_sca])

    for j in trange(Nmc+1,
                    desc='IRS SL and LL matching [MC]'):
        for i, ch in enumerate(chnl[:4]):
            if j==0:
                ic = intercalib(path_tmp+src+'_'+ch+'_conv')
                ic.correct_spec(filOUT=path_tmp+src+'_'+ch+'_match',
                                gain=scale[i])
            else:
                ic = intercalib(path_tmp+src+'_'+ch+'_'+str(j)+'_conv')
                ic.correct_spec(filOUT=path_tmp+src+'_'+ch+'_'+str(j)+'_match',
                                gain=scale[i])
    
## Match SH-LH
##-------------
if hires_match=='y':
    ## SH
    data_sh = read_fits(path_tmp+src+'_SH_conv').data
    wvl_sh = read_fits(path_tmp+src+'_SH_conv').wave
    ## LH
    data_lh = read_fits(path_tmp+src+'_LH_conv').data
    wvl_lh = read_fits(path_tmp+src+'_LH_conv').wave
    
    iwmax_sh = closest(wvl_lh, wvl_sh[-1], side='left') + 1 # LH index
    iwmin_lh = closest(wvl_sh, wvl_lh[0], side='right') # SH index
    right_sh = trapz(data_sh[iwmin_lh:], wvl_sh[iwmin_lh:],
                     dx=wvl_sh[1]-wvl_sh[0], axis=0)
    left_lh = trapz(data_lh[:iwmax_sh], wvl_lh[:iwmax_sh],
                    dx=wvl_lh[1]-wvl_lh[0], axis=0)

    ## Integrated scaling factor
    scale = right_sh / left_lh
    
    ## One point scaling factor
    # i_sh = closest(wvl_sh, 19.10)
    # i_lh = closest(wvl_lh, 19.30)
    # scale = data_sh[i_sh] / data_lh[i_lh]
    
    ## Display scale map
    mask_sca = ~np.isnan(scale)
    print('LH to SH scaling factor: ', scale[mask_sca])

    for j in trange(Nmc+1,
                    desc='IRS SH-LH matching [MC]'):
        if j==0:
            ic = intercalib(path_tmp+src+'_LH'+'_conv')
            ic.correct_spec(filOUT=path_tmp+src+'_LH'+'_match',
                            gain=scale)
        else:
            ic = intercalib(path_tmp+src+'_LH'+'_'+str(j)+'_conv')
            ic.correct_spec(filOUT=path_tmp+src+'_LH'+'_'+str(j)+'_match',
                            gain=scale)
                
if concat_irs=='y':

    ## Define wrange for SH and LH
    if ('SH' in chnl) or ('LH' in chnl):
        wrange = [ (10.00, 19.19), # sh
                   (19.20, 37.10), ] # lh
    else:
        wrange = None
    
    for j in trange(Nmc+1,
                    desc='IRS concatenation [MC]'):
        files = []
        for i in range(Nch_s):
            if j==0:
                if lores_match=='y':
                    files.append(path_tmp+src+'_'+chnl[i]+'_match')
                elif hires_match=='y':
                    print('coucou')
                    if chnl[i]=='LH':
                        files.append(path_tmp+src+'_'+chnl[i]+'_match')
                    else:
                        files.append(path_tmp+src+'_'+chnl[i]+'_conv')
                else:
                    files.append(path_tmp+src+'_'+chnl[i]+'_conv')
            else:
                if lores_match=='y':
                    files.append(path_tmp+src+'_'+chnl[i]+'_'+str(j)+'_match')
                elif hires_match=='y':
                    if chnl[i]=='LH':
                        files.append(path_tmp+src+'_'+chnl[i]+'_'+str(j)+'_match')
                    else:
                        files.append(path_tmp+src+'_'+chnl[i]+'_'+str(j)+'_conv')
                else:
                    files.append(path_tmp+src+'_'+chnl[i]+'_'+str(j)+'_conv')

        concatenate(files, out_irs+'_'+str(j), wsort=False, wrange=wrange,
                    keepfrag=keepfrag, cropedge=cropedge)
    
    ## Cal unc (all)
    ##---------------
    mcimage = []
    for j in trange(Nmc+1, leave=False,
                    desc='IRS Reading [MC]'):
        if j==0:
            hd = read_fits(out_irs+'_0')
            header = hd.header
            wvl = hd.wave
        else:
            hd = read_fits(out_irs+'_'+str(j))
            mcimage.append(hd.data)
    if Nmc>1:
        mcimage = np.array(mcimage)
        unc = np.nanstd(mcimage, axis=0)
        write_fits(out_irs+'_unc', header, unc, wvl)

##----------------------------------------------------------

##                   Inter-Calibration

##----------------------------------------------------------
pre_calib1 = input("Run inter-calibration prepipeline (DustPedia) (y/n)? ")
int_calib = input("Run inter-calibration (y/n)? ")

if not os.path.exists(path_tmp+'calib/'):
    os.makedirs(path_tmp+'calib/')
    
for i in range(Nphot):

    ## Prepare photometry
    ##--------------------
    if phot[i]=='IRAC4':
        phot_ker = path_ker+'Kernel_HiRes_IRAC_8.0_to_Gauss_06.0'
    elif phot[i]=='MIPS1':
        phot_ker = path_ker+'Kernel_HiRes_MIPS_24_to_Gauss_06.0'
    
    raw_phot1 = path_phot+src+'_'+phot[i]+'_DP'
    fits_phot1 = path_cal+src+'_'+phot[i]+'_DP'
    tmp_phot1 = path_tmp+'calib/'+src+'_'+phot[i]+'_DP'
    
    fits_spec = path_cal+src+'_'+phot[i]+'_IRS'
    tmp_spec = path_tmp+'calib/'+src+'_'+phot[i]+'_IRS'
    
    ## DustPedia (phot1)
    ##===================
    if pre_calib1=='y':

        ## Convert phot unit (Jy/pix -> MJy/sr)
        ## This step should be before reprojection (pixscale changed)
        Jy_per_pix_to_MJy_per_sr(raw_phot1, filOUT=path_tmp+src+'_'+phot[i]+'_DP')
        Jy_per_pix_to_MJy_per_sr(raw_phot1+'_unc',
                                 filOUT=path_tmp+src+'_'+phot[i]+'_DP_unc')
        
        ## Reproject phot
        coadd_fp_crop = fixwcs(out_irs+'_0'+fitsext).header # if cropedge
        mtg = imontage('interp', tmpdir=path_tmp) # In 'exact' case WCS failed to converge
        mtg.reproject_mc(path_tmp+src+'_'+phot[i]+'_DP', refheader=coadd_fp_crop,
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

    ## Fit correlation
    ##-----------------
    ## Data
    data_phot1 = read_fits(fits_phot1).data
    unc_phot1 = read_fits(fits_phot1+'_unc').data
    data_spec = read_fits(fits_spec).data
    unc_spec = read_fits(fits_spec+'_unc').data
    if i==0:
        d1_phot1 = data_phot1
        u1_phot1 = unc_phot1
        d1_spec = data_spec
        u1_spec = unc_spec
    elif i==1:
        d2_phot1 = data_phot1
        u2_phot1 = unc_phot1
        d2_spec = data_spec
        u2_spec = unc_spec
            
    if int_calib=='y':
        
        ## Plot
        pix_phot1 = data_phot1[:,:].reshape((-1,))
        pix_spec = data_spec[:,:].reshape((-1,))
        pix_phot1_unc = unc_phot1[:,:].reshape((-1,))
        pix_spec_unc = unc_spec[:,:].reshape((-1,))
        
        # xgrid = np.arange(0,1e3,1)
        xgrid = np.logspace(-6,4,10000)
        
        ## NaN mask
        mask = ~np.logical_or( np.isnan(pix_spec), np.isnan(pix_phot1) )
    
        # if i==0:
        #     mask_all = mask
        # else:
        #     mask_all = np.logical_and( mask, mask_all )
        
        ## S/N ratio
        # print('S/N ('+phot[i]+' - IRS) = \n', pix_spec[mask]/pix_spec_unc[mask])
        # exit()
            
        ## IRS - DustPedia
        p0 = pplot(pix_phot1, pix_spec,
                   yerr=pix_spec_unc, xerr=pix_phot1_unc,
                   fmt='.', c='k', ec='r', legend='upper left',
                   xlog=1, ylog=1, nonposx='clip', nonposy='clip',
                   xlim=(1e-1,1e4), ylim=(1e-1,1e4),
                   xlab='DustPedia (MJy/sr)', ylab='IRS (MJy/sr)',
                   figsize=(8,8), title=src+'_'+phot[i])
        
        ## Linear fit
        if i==0:
            popt, pcov = curve_fit(f_lin0, pix_spec[mask], pix_phot1[mask])
                                   # sigma=pix_phot1_unc[mask])
            calib_gain = popt[0]
            # calib_off = popt[1]
            calib_off = 0.
            # calib_gain = 1.
            # calib_off = popt[0]
            print('Inter-Calibration ('+phot[i]+') gain = {:.4}'.format(calib_gain))
            print('Inter-Calibration ('+phot[i]+') offset = {:.4}'.format(calib_off))
            label = phot[i]+': y={0:.4}x+{1:.4}'.format(calib_gain, calib_off)
            p.add_plot(f_lin0(xgrid, *popt), xgrid,
                       c='y', ls='-', label=label)
            
            # popt, pcov = curve_fit(f_lin0, pix_phot1[mask], pix_spec[mask],)
            #                        # sigma=pix_spec_unc[mask])
            # calib_gain = 1. / popt[0]
            # # calib_off = -popt[1] / popt[0]
            # calib_off = 0.
            # print('Inter-Calibration ('+phot[i]+') gain = {:.4}'.format(calib_gain))
            # print('Inter-Calibration ('+phot[i]+') offset = {:.4}'.format(calib_off))
            # label = phot[i]+': y={0:.4}x+{1:.4}'.format(calib_gain, calib_off)
            # p.add_plot(xgrid, f_lin0(xgrid, *popt),
            #            c='y', ls='-', label=label)
        else:
            # mask_mips1 = ~pix_phot1<1.e0
            # mask = np.logical_and( mask, mask_mips1 )
            
            popt, pcov = curve_fit(f_lin0, pix_spec[mask], pix_phot1[mask])
                                   # sigma=pix_phot1_unc[mask])
            calib_gain = popt[0]
            # calib_off = popt[1]
            calib_off = 0.
            print('Inter-Calibration ('+phot[i]+') gain = {:.4}'.format(calib_gain))
            print('Inter-Calibration ('+phot[i]+') offset = {:.4}'.format(calib_off))
            label = phot[i]+': y={0:.4}x+{1:.4}'.format(calib_gain, calib_off)
            p.add_plot(f_lin0(xgrid, *popt), xgrid,
                       c='y', ls='-', label=label)

            # popt, pcov = curve_fit(f_lin0, pix_phot1[mask], pix_spec[mask],)
            #                        # sigma=pix_spec_unc[mask])
            # calib_gain = 1. / popt[0]
            # # calib_off = -popt[1] / popt[0]
            # calib_off = 0.
            # print('Inter-Calibration ('+phot[i]+') gain = {:.4}'.format(calib_gain))
            # print('Inter-Calibration ('+phot[i]+') offset = {:.4}'.format(calib_off))
            # label = phot[i]+': y={0:.4}x+{1:.4}'.format(calib_gain, calib_off)
            # p.add_plot(xgrid, f_lin0(xgrid, *popt),
            #            c='y', ls='-', label=label)
            
        p.set_font()
        
        p.save(path_cal+'intercalib_'+phot[i])
        
        
        ## Spectral correction
        ##---------------------
        sc = intercalib(out_irs+'_0')
        if i==0:
            sc.correct_spec(filOUT=out_irs) # Duplicate spectral cube
        sc.read_filter(phot[i]) # Convert calib_off to spec_off
        spec_off = calib_off * sc.specoff_ov_bboff[0]
        spec_gain = calib_gain
        if ('SH' in chnl) or ('LH' in chnl):
            wmin = 9.0
            wmax = 39.00
        else:
            if i==0:
                wmin = 5.13
                wmax = 14.28
            elif i==1:
                wmin = 14.29
                wmax = 39.00
            
        do_correct = input(" - Do IRS spectral correction [{}] (y/n)? ".format(phot[i]))
        if do_correct=='y':
            do_indiv = input("   - Do individual correction for each pixel (y/n)? ")

            if do_indiv=='y':
                ## Individual
                ##============
                spec_gain = (data_phot1-spec_off)/data_spec
                mask2D = ~np.logical_or( np.isnan(data_spec), np.isnan(data_phot1) )
                imp = improve(out_irs)
                im = imp.im
                for k, lam in enumerate(imp.wvl):
                    if lam>=wmin and lam<=wmax:
                        im[k][mask2D] = imp.im[k][mask2D] * spec_gain[mask2D] + spec_off
                ## Supplementary mask
                # mask2D_s = (data_spec/data_phot1>1.) # Define supp creterion
                # imp0 = improve(out_irs+'_0')
                # for k, lam in enumerate(imp.wvl):
                #     im[k][mask2D_s] = imp0.im[k][mask2D_s]
                
                write_fits(out_irs, imp.hdr, im, imp.wvl, imp.wmod)
                
            else:
                ## Global
                ##========
                sc = intercalib(out_irs)
                sc.correct_spec(filOUT=out_irs,
                                gain=spec_gain, offset=spec_off,
                                wlim=(wmin,wmax))

##----------------------------------------------------------

##                      Plot spectra

##----------------------------------------------------------
do_plot = input("Plot IRS spectra (y/n)? ")
if do_plot=='y':
    ds = read_fits(out_irs, out_irs+'_unc')
    # ds.data = respect().smooth(out_irs, out_irs, lim_unc=1.e2, cmin=2)
    ds0 = read_fits(out_irs+'_0')
    # ds_irc = read_fits(out_irc)
    pp = intercalib(out_irs)
    pp.read_filter(phot)
    wcen = pp.wcen
    Nw, Ny, Nx = ds.data.shape
    for x in range(Nx):
        for y in range(Ny):
            if ~np.isnan(ds.data[:,y,x]).any():
                # if ~np.isnan(ds_irc.data[0,y,x]).any():
                p = pplot(ds.wave, ds.data[:,y,x], yerr=ds.unc[:,y,x],
                          c='k', lw=.7, ec='r', xlog=1, ylog=1, legend='upper left')                
                p.add_plot(ds0.wave, ds0.data[:,y,x],
                           c='y', ls='--', zorder=-1)
                
                if len(phot)>0:
                    p.add_plot(wcen[0], d1_phot1[y,x], yerr=u1_phot1[y,x],
                               c='m', marker='o', ms=10, zorder=100, label=phot[0])
                    p.add_plot(wcen[0], d1_spec[y,x], yerr=u1_spec[y,x],
                               c='g', marker='^', ms=10, zorder=101, label='IRS-'+phot[0])
                if len(phot)>1:
                    p.add_plot(wcen[1], d2_phot1[y,x], yerr=u1_phot1[y,x],
                               c='orange', marker='o', ms=10, zorder=102, label=phot[1])
                    p.add_plot(wcen[1], d2_spec[y,x], yerr=u1_spec[y,x],
                               c='c', marker='^', ms=10, zorder=103, label='IRS-'+phot[1])

                if ('SH' in chnl) or ('LH' in chnl):
                    xtic = [10, 12, 15, 20, 30, 40]
                    xtic_min = np.arange(9., 41., 2.)
                    ## Zoom
                    # xtic = np.arange(19., 19.4, .1)
                    # xtic_min = np.arange(19., 19.4, .02)
                    # p.set_ax(xlim=(19,19.4))
                else:
                    xtic = [5, 6, 7, 8, 10, 12, 15, 20, 30, 40]
                    xtic_min = np.arange(5., 41., 2.)
                p.ax.set_xticks(xtic, minor=False) # major
                p.ax.set_xticks(xtic_min, minor=True) # minor
                p.ax.xaxis.set_major_formatter(ScalarFormatter()) # major
                p.ax.xaxis.set_minor_formatter(NullFormatter()) # minor
                
                p.save(path_fig+'IRS_('+str(x+1)+','+str(y+1)+')')
