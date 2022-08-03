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
from scipy.optimize import curve_fit
from matplotlib.ticker import ScalarFormatter, NullFormatter

## laputan
from laputan.inout import fclean, fitsext, read_fits, write_fits
from laputan.imaging import ( concatenate, iconvolve, iswarp, iuncert,
                              improve, respect, Jy_per_pix_to_MJy_per_sr )
from laputan.astrom import fixwcs
from laputan.calib import intercalib
from laputan.maths import f_lin, f_lin0
from laputan.plots import pplot

##---------------------------
##      Initialisation
##---------------------------
## Local
from buildinfo import ( src, Nmc, verbose,
                        path_phot, path_ker, path_cal,
                        path_idl, csv_ker, path_tmp,
                        path_out, out_irc, out_irs, path_fig
)
# Nmc = 2

##---------------------------
##   Concatenate IRC & IRS
##---------------------------
concatenate([out_irc,out_irs], path_tmp+src+'_ma', wsort=False)
data0 = read_fits(path_tmp+src+'_ma').data
ma = np.ma.array(data0, mask=np.isnan(data0))
mask_any = ma.mask.any(axis=0)
mask_all = ma.mask.all(axis=0)

concatenate([out_irc,out_irs], path_out+src, wsort=False)
concatenate([out_irc+'_unc',out_irs+'_unc'], path_out+src+'_unc', wsort=False)
hd = read_fits(path_out+src, path_out+src+'_unc')
data = hd.data
unc = hd.unc
header = hd.header
wvl = hd.wave
for k in range(len(wvl)):
    data[k][mask_all] = np.nan
    unc[k][mask_all] = np.nan
write_fits(path_out+src, header, data, wvl)
write_fits(path_out+src+'_unc', header, unc, wvl)

##---------------------------
##     Inter-Calibration
##---------------------------
## Reprojection footprint
refheader = fixwcs(path_out+src+fitsext).header
swp = iswarp(refheader=refheader,
             tmpdir=path_tmp, verbose=verbose)

Nx = refheader['NAXIS1']
Ny = refheader['NAXIS2']

phot = ['IRAC2', 'IRAC3']
Nphot = len(phot)

## Prepare photometry
##--------------------
precalib = input("Run inter-calibration prepipeline (y/n)? ")

for i in range(Nphot):
    
    if phot[i]=='IRAC2':
        phot_ker = path_ker+'Kernel_HiRes_IRAC_4.5_to_Gauss_06.0'
    elif phot[i]=='IRAC3':
        phot_ker = path_ker+'Kernel_HiRes_IRAC_5.8_to_Gauss_06.0'
    
    raw_phot1 = path_phot+src+'_'+phot[i]+'_DP'
    fits_phot1 = path_cal+src+'_'+phot[i]+'_DP'
    
    raw_phot2 = path_phot+src+'_'+phot[i]+'_SINGS'
    wgt_phot2 = path_phot+src+'_'+phot[i]+'_SINGS_wgt'
    fits_phot2 = path_cal+src+'_'+phot[i]+'_SINGS'
    
    fits_spec = path_cal+src+'_'+phot[i]+'_MIR'
    
    if precalib:
    
        ## DustPedia (phot1)
        ##===================
        
        ## Convert phot unit (Jy/pix -> MJy/sr)
        ## This step should be before iswarp reprojection (pixscale changed)
        Jy_per_pix_to_MJy_per_sr(raw_phot1, filOUT=fits_phot1)
        if i==1:
            Jy_per_pix_to_MJy_per_sr(raw_phot1+'_unc', filOUT=fits_phot1+'_unc')

            for j in trange(Nmc+1, leave=False,
                    desc='Preparing photometry [MC]'):
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
            
            ## Propagate DustPedia uncertainty
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
        else:
            ## Reproject phot
            comb = swp.combine(fits_phot1,# combtype='wgt_avg',
                               cropedge=False, tmpdir=path_tmp,
                               filOUT=fits_phot1).image

            ## Convolve phot
            conv = iconvolve(fits_phot1, phot_ker, csv_ker, filOUT=fits_phot1)
            
            conv.do_conv(idldir=path_idl)
        
        ## SINGS (phot2)
        ##===============
        
        ## Create uncertainty via weight map (suppose uniform contribution)
        ## The weight maps contain the information on the number of frames 
        ## that were used to create the science mosaics at each pixel (value= # frames x10); 
        ## the pixel size of the weight maps is the same as the science mosaics.
        ## -- SINGS v5 release doc
        if phot[i]=='IRAC2':
            bg = read_fits(raw_phot2).data[542:582,1102:1142]
            bg_wgt = read_fits(wgt_phot2).data[542:582,1102:1142] / 10.
        elif phot[i]=='IRAC3':
            bg = read_fits(raw_phot2).data[402:442,1022:1062]
            bg_wgt = read_fits(wgt_phot2).data[402:442,1022:1062] / 10.
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
        for j in trange(Nmc+1,
                        desc='MIR sythetic photometry [MC]'):
            if j==0:
                ic0 = intercalib(path_out+src)
                sp0 = ic0.synthetic_photometry(phot[i])
                write_fits(fits_spec, ic0.hdr, sp0.Fnu_filt)
            else:
                imp = improve(path_out+src)
                imp.rand_norm(path_out+src+'_unc')
                write_fits(path_tmp+src+'_'+str(j), imp.hdr, imp.im, imp.wvl)
                ic = intercalib(path_tmp+src+'_'+str(j))
                sp = ic.synthetic_photometry(phot[i])
                ## Cal unc
                mcimage.append(sp.Fnu_filt)
        if Nmc>1:
            print('Calculating uncertainty cube...')
            mcimage = np.array(mcimage)
            unc = np.nanstd(mcimage, axis=0)
            write_fits(fits_spec+'_unc', ic0.hdr, unc)

    ## Inter-Calibration
    ##-------------------
    data_phot1 = read_fits(fits_phot1).data
    pix_phot1 = data_phot1[:,:].reshape((-1,))
    # if i==1:
    #     unc_phot1 = read_fits(fits_phot1+'_unc').data
        
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
    # print('S/N ('+phot[i]+' - MIR) = \n', pix_spec[mask]/pix_spec_unc[mask])
    # print('S/N ('+phot[i]+' - SINGS) = \n', pix_phot2[mask]/pix_phot2_unc[mask])
    # exit()
        
    ## SINGS - DP
    p0 = pplot(pix_phot1,  pix_phot2, yerr=pix_phot2_unc,
               fmt='.', c='y', ec='r',elw=10,
               xlog=1, ylog=1,
               xlim=(0,1e3), ylim=(0,1e3),
               xlab='DustPedia (MJy/sr)', ylab='SINGS (MJy/sr)',
               figsize=(8,8), title='M82 '+phot[i])
    p0.save(path_cal+'SINGS-DP_'+phot[i])
        
    ## MIR spec - SINGS
    p = pplot(fmt='.', legend='upper left',
              xlog=1, ylog=1,
              xlim=(0,1e3), ylim=(0,1e3),
              xlab='MIR spec (MJy/sr)', ylab='SINGS (MJy/sr)', 
              figsize=(8,8), title='M82 '+phot[i]+' inter-calibration')
    ## MIR spec - SINGS
    p.add_plot(pix_spec, pix_phot2, xerr=pix_spec_unc, yerr=pix_phot2_unc,
               fmt='.', c='y', ec='r')
    
    ## Linear fit
    popt, pcov = curve_fit(f_lin, pix_spec[mask], pix_phot2[mask],
                           sigma=pix_phot2_unc[mask])
    calib_fac = popt[0]
    calib_off = popt[1]
    # popt, pcov = curve_fit(f_lin0, pix_phot2[mask], pix_spec[mask])
    # calib_fac = 1. / popt[0]
    # calib_off = -popt[1] / popt[0]
    
    print('Inter-Calibration ('+phot[i]+') factor = {:.4}'.format(calib_fac))
    print('Inter-Calibration ('+phot[i]+') offset = {:.4}'.format(calib_off))
    # label = src+' calib fac = {:.4}'.format(calib_fac)
    label = src+': y={0:.4}x+{1:.4}'.format(calib_fac, calib_off)
    
    p.add_plot(xgrid, f_lin(xgrid, *popt),
               c='k', ls='-', label=label)
    
    p.set_font()
    
    p.save(path_cal+'calib_'+phot[i])

##---------------------------
##       Plot spectra
##---------------------------
do_plot = input("Plot MIR spectra (y/n)? ")
if do_plot:
    ds1 = read_fits(out_irc+'_0')
    ds2 = read_fits(out_irs+'_0')
    ds = read_fits(path_out+src, path_out+src+'_unc')
    # ds.data = respect().smooth(path_out+src, path_out+src, lim_unc=1.e2, cmin=2)
    wcen = intercalib().wcenter(phot)
    for x in range(Nx):
        for y in range(Ny):
            if ~np.isnan(ds.data[:,y,x]).any():
                p = pplot(ds.wave, ds.data[:,y,x], yerr=ds.unc[:,y,x],
                          c='k', lw=.7, ec='r', xlog=1, ylog=1,
                          legend='upper left', title=src+'_('+str(x)+','+str(y)+')')
                p.add_plot(wcen[0], d1_phot2[y,x], yerr=u1_phot2[y,x],
                           c='m', marker='o', ms=10, zorder=100, label=phot[0])
                p.add_plot(wcen[0], d1_spec[y,x], yerr=u1_spec[y,x],
                           c='g', marker='^', ms=10, zorder=101, label='MIR spec-'+phot[0])
                p.add_plot(wcen[1], d2_phot2[y,x], yerr=u1_phot2[y,x],
                           c='orange', marker='o', ms=10, zorder=102, label=phot[1])
                p.add_plot(wcen[1], d2_spec[y,x], yerr=u1_spec[y,x],
                           c='c', marker='^', ms=10, zorder=103, label='MIR spec-'+phot[1])


                ph = ['IRAC1','IRAC4','MIPS1']
                wc = intercalib().wcenter(ph)
                c1 = ['pink','tomato','maroon']
                c2 = ['lightgreen','lime','olive']
                for i,ip in enumerate(ph):
                    data_phot = read_fits(path_cal+src+'_'+ip+'_SINGS').data
                    unc_phot = read_fits(path_cal+src+'_'+ip+'_SINGS_unc').data
                    if i==0:
                        data_spec = read_fits(path_cal+src+'_'+ip+'_IRC').data
                        unc_spec = read_fits(path_cal+src+'_'+ip+'_IRC_unc').data
                    else:
                        data_spec = read_fits(path_cal+src+'_'+ip+'_IRS').data
                        unc_spec = read_fits(path_cal+src+'_'+ip+'_IRS_unc').data
                    p.add_plot(wc[i], data_phot2[y,x], yerr=unc_phot2[y,x],
                               c=c1[i], marker='o', ms=10, zorder=102, label=ph[i])
                    if i==0:
                        p.add_plot(wc[i], data_spec[y,x], yerr=unc_spec[y,x],
                                   c=c2[i], marker='^', ms=10, zorder=103, label='IRC-'+ph[i])
                    else:
                        p.add_plot(wc[i], data_spec[y,x], yerr=unc_spec[y,x],
                                   c=c2[i], marker='^', ms=10, zorder=103, label='IRS-'+ph[i])
                
                p.add_plot(ds1.wave, ds1.data[:,y,x],
                           c='y', lw=.7, ls='-', zorder=-1, label='IRC raw')
                p.add_plot(ds2.wave, ds2.data[:,y,x],
                           c='y', lw=.7, ls='-', zorder=-1, label='IRS raw')

                # xtic = [2, 3, 4, 5, 6, 7, 8, 10, 12, 15, 20, 30, 40]
                xtic = []
                xtic_min = np.arange(2., 41., 1.)
                p.ax.set_xticks(xtic, minor=False) # major
                p.ax.set_xticks(xtic_min, minor=True) # minor
                p.ax.xaxis.set_major_formatter(ScalarFormatter()) # major
                p.ax.xaxis.set_minor_formatter(NullFormatter()) # minor
                
                p.save(path_fig+src+'_('+str(x)+','+str(y)+')', transparent=1)
                # p.save(path_out+src+'_(60,22)')
