#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""

Build stitched AKARI/IRC-Spitzer/IRS spectral cube

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
from buildinfo import ( src, Nmc, verbose,
                        path_idl, path_phot, path_cal, path_ker, csv_ker, 
                        out_irc, out_irs, path_tmp, path_out, path_fig
)
Nmc = 20

phot = ['IRAC2', 'IRAC3']
Nphot = len(phot)

##---------------------------
##   Concatenate IRC & IRS
##---------------------------
concat_mir = input("Stitch MIR spectra (y/n)? ")
if concat_mir=='y':
    concatenate([out_irc,out_irs], path_tmp+src+'_ma', wsort=False)
    data0 = read_fits(path_tmp+src+'_ma').data
    ma = np.ma.array(data0, mask=np.isnan(data0))
    mask_any = ma.mask.any(axis=0)
    mask_all = ma.mask.all(axis=0)
    
    Nw, Ny, Nx = data0.shape
    
    for j in trange(Nmc+1,
                    desc='MIR concatenation [MC]'):
        if j==0:
            concatenate([out_irc,out_irs], path_out+src+'_0', wsort=False)
            concatenate([out_irc+'_unc',out_irs+'_unc'], path_out+src+'_unc', wsort=False)
        else:
            concatenate([out_irc+'_'+str(j),out_irs+'_'+str(j)],
                        path_out+src+'_'+str(j), wsort=False)
        hd = read_fits(path_out+src+'_'+str(j))
        data = hd.data
        header = hd.header
        wvl = hd.wave
        for k in range(Nw):
            data[k][mask_all] = np.nan
        write_fits(path_out+src+'_'+str(j), header, data, wvl)
        if j==0:
            unc = read_fits(path_out+src+'_unc').data
            for k in range(Nw):
                unc[k][mask_all] = np.nan
            write_fits(path_out+src+'_unc', header, unc, wvl)

##---------------------------
##    Prepare photometry
##---------------------------
pre_calib1 = input("Run inter-calibration prepipeline (DustPedia) (y/n)? ")
pre_calib2 = input("Run inter-calibration prepipeline (SINGS) (y/n)? ")
int_calib = input("Run inter-calibration (y/n)? ")

if not os.path.exists(path_tmp+'calib/'):
    os.makedirs(path_tmp+'calib/')
    
for i in range(Nphot):
    
    if phot[i]=='IRAC2':
        phot_ker = path_ker+'Kernel_HiRes_IRAC_4.5_to_Gauss_06.0'
    elif phot[i]=='IRAC3':
        phot_ker = path_ker+'Kernel_HiRes_IRAC_5.8_to_Gauss_06.0'
    
    raw_phot1 = path_phot+src+'_'+phot[i]+'_DP'
    fits_phot1 = path_cal+src+'_'+phot[i]+'_DP'
    tmp_phot1 = path_tmp+'calib/'+src+'_'+phot[i]+'_DP'
    
    raw_phot2 = path_phot+src+'_'+phot[i]+'_SINGS'
    wgt_phot2 = path_phot+src+'_'+phot[i]+'_SINGS_wgt'
    fits_phot2 = path_cal+src+'_'+phot[i]+'_SINGS'
    tmp_phot2 = path_tmp+'calib/'+src+'_'+phot[i]+'_SINGS'
    
    fits_spec = path_cal+src+'_'+phot[i]+'_MIR'
    tmp_spec = path_tmp+'calib/'+src+'_'+phot[i]+'_MIR'
    
    ## DustPedia (phot1)
    ##===================
    if pre_calib1=='y':
        
        ## Convert phot unit (Jy/pix -> MJy/sr)
        ## This step should be before iswarp reprojection (pixscale changed)
        Jy_per_pix_to_MJy_per_sr(raw_phot1, filOUT=path_tmp+src+'_'+phot[i]+'_DP')

        if i==0:
            ## DustPedia IRAC2 of M82 has no uncertainty map
            ## Create homogeneous uncertainty map
            bg = read_fits(path_tmp+src+'_'+phot[i]+'_DP').data[1720:1760,3110:3150]
            iuncert(path_tmp+src+'_'+phot[i]+'_DP',
                    filOUT=path_tmp+src+'_'+phot[i]+'_DP_unc', BG_image=bg)
        else:
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
        # if phot[i]=='IRAC2':
        #     bg = read_fits(raw_phot2).data[542:582,1102:1142]
        #     bg_wgt = read_fits(wgt_phot2).data[542:582,1102:1142] / 10.
        # elif phot[i]=='IRAC3':
        #     bg = read_fits(raw_phot2).data[402:442,1022:1062]
        #     bg_wgt = read_fits(wgt_phot2).data[402:442,1022:1062] / 10.
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
        for j in trange(Nmc+1, leave=False,
                        desc='MIR sythetic photometry [MC]'):
            ic = intercalib(path_out+src+'_'+str(j))
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
        # print('S/N ('+phot[i]+' - MIR) = \n', pix_spec[mask]/pix_spec_unc[mask])
        # print('S/N ('+phot[i]+' - SINGS) = \n', pix_phot2[mask]/pix_phot2_unc[mask])
        # exit()
            
        ## SINGS - DP
        p0 = pplot(pix_phot1, pix_phot2,
                   yerr=pix_phot2_unc, xerr=pix_phot1_unc,
                   fmt='.', c='k', ec='r',elw=1,
                   xlog=1, ylog=1,
                   xlim=(0,1e4), ylim=(0,1e4),
                   xlab='DustPedia (MJy/sr)', ylab='SINGS (MJy/sr)',
                   figsize=(8,8), title='M82 '+phot[i])
        p0.save(path_cal+'SINGS-DP_'+phot[i])
            
        ## MIR spec - SINGS
        p = pplot(pix_phot2, pix_spec,
                  yerr=pix_spec_unc, xerr=pix_phot2_unc,
                  fmt='.', c='k', ec='r', legend='upper left',
                  xlog=1, ylog=1, nonposx='clip', nonposy='clip',
                  xlim=(1e-6,1e4), ylim=(1e-6,1e4),
                  xlab='SINGS (MJy/sr)', ylab='MIR spec (MJy/sr)',
                  figsize=(8,8), title=src+'_'+phot[i]+' inter-calibration')
        
        ## Linear fit
        # popt, pcov = curve_fit(f_lin, pix_spec[mask], pix_phot2[mask],
        #                        sigma=pix_phot2_unc[mask])
        # calib_fac = popt[0]
        # calib_off = popt[1]

        if i==1:
            mask_irac3 = ~np.ma.masked_where(pix_spec/pix_phot2<1.e-2, pix_phot2).mask
            mask = np.ma.array(mask,
                               mask=np.logical_and(mask,mask_irac3)).mask
            
        popt, pcov = curve_fit(f_lin0, pix_phot2[mask], pix_spec[mask],
                               sigma=pix_spec_unc[mask])
        calib_fac = 1. / popt[0]
        # calib_off = -popt[1] / popt[0]
        
        print('Inter-Calibration ('+phot[i]+') factor = {:.4}'.format(calib_fac))
        # print('Inter-Calibration ('+phot[i]+') offset = {:.4}'.format(calib_off))
        label = phot[i]+' calib fac = {:.4}'.format(calib_fac)
        # label = phot[i]+': y={0:.4}x+{1:.4}'.format(calib_fac, calib_off)
        
        p.add_plot(xgrid, f_lin0(xgrid, *popt),
                   c='y', ls='-', label=label)
        
        p.set_font()
        
        p.save(path_cal+'intercalib_'+phot[i])

##---------------------------
##       Plot spectra
##---------------------------
do_plot = input("Plot MIR spectra (y/n)? ")
if do_plot:
    ds1 = read_fits(out_irc+'_0')
    ds2 = read_fits(out_irs+'_0')
    ds = read_fits(path_out+src+'_0', path_out+src+'_unc')
    Nw, Ny, Nx = ds.data.shape
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

                xtic = [2, 3, 4, 5, 6, 7, 8, 10, 12, 15, 20, 30, 40]
                # xtic = []
                xtic_min = np.arange(2., 41., 1.)
                p.ax.set_xticks(xtic, minor=False) # major
                p.ax.set_xticks(xtic_min, minor=True) # minor
                p.ax.xaxis.set_major_formatter(ScalarFormatter()) # major
                p.ax.xaxis.set_minor_formatter(NullFormatter()) # minor
                
                p.save(path_fig+src+'_('+str(x)+','+str(y)+')', transparent=1)
