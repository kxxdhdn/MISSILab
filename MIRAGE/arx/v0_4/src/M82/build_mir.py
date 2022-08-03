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

## rapyuta
from rapyuta.inout import fclean, fitsext, read_fits, write_fits
from rapyuta.imaging import ( iconvolve, iswarp, iuncert, concatenate, icrop,
                              imontage, improve, Jy_per_pix_to_MJy_per_sr )
from rapyuta.astrom import fixwcs
from rapyuta.arrays import closest
from rapyuta.calib import intercalib
from rapyuta.maths import f_lin, f_lin0, f_lin1
from rapyuta.plots import pplot

##----------------------------------------------------------

##                     Preliminaries

##----------------------------------------------------------

## Local
from buildinfo import ( src, Nmc, verbose, coadd_fp,
                        path_idl, path_phot, path_cal, path_ker, csv_ker, 
                        out_irc, out_irs, path_tmp, path_out, path_fig
)

## Banner
print('\n============================================================\n')

print('        MIRAGE - Mid-infrared cube builder - '+src)

print('\n============================================================\n')

# Nmc = 2

phot = ['IRAC2', 'IRAC3']
Nphot = len(phot)

##----------------------------------------------------------

##                 Concatenate IRC & IRS

##----------------------------------------------------------
concat_mir = input("Stitch MIR spectra (y/n)? ")
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
    
if concat_mir=='y':
    concatenate([out_irc,out_irs], path_tmp+src+'_ma', wsort=False)
    data0 = read_fits(path_tmp+src+'_ma').data
    mask_any = np.isnan(data0).any(axis=0)
    mask_all = np.isnan(data0).all(axis=0)
    
    Nw, Ny, Nx = data0.shape
    for j in trange(Nmc+1,
                    desc='MIR concatenation [MC]'):
        if j==0:
            concatenate([out_irc,out_irs],
                        path_out+src+'_0', wsort=False,
                        keepfrag=keepfrag, cropedge=cropedge)
            concatenate([out_irc+'_unc',out_irs+'_unc'],
                        path_out+src+'_unc', wsort=False,
                        keepfrag=keepfrag, cropedge=cropedge)
        else:
            concatenate([out_irc+'_'+str(j),out_irs+'_'+str(j)],
                        path_out+src+'_'+str(j), wsort=False,
                        keepfrag=keepfrag, cropedge=cropedge)
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

##----------------------------------------------------------

##                   Inter-Calibration

##----------------------------------------------------------
pre_calib1 = input("Run inter-calibration prepipeline (DustPedia) (y/n)? ")
pre_calib2 = input("Run inter-calibration prepipeline (SINGS) (y/n)? ")
int_calib = input("Run inter-calibration (y/n)? ")

if not os.path.exists(path_tmp+'calib/'):
    os.makedirs(path_tmp+'calib/')
    
for i in range(Nphot):

    ## Prepare photometry
    ##--------------------
    if phot[i]=='IRAC2':
        phot_ker = path_ker+'Kernel_HiRes_IRAC_4.5_to_Gauss_06.0'
    elif phot[i]=='IRAC3':
        phot_ker = path_ker+'Kernel_HiRes_IRAC_5.8_to_Gauss_06.0'
    
    raw_phot1 = path_phot+src+'_'+phot[i]+'_DP'
    fits_phot1 = path_cal+src+'_'+phot[i]+'_DP'
    tmp_phot1 = path_tmp+'calib/'+src+'_'+phot[i]+'_DP'
    
    raw_phot2 = path_phot+src+'_'+phot[i]+'_SINGS'
    wgt_phot2 = path_phot+src+'_'+phot[i]+'_SINGS_wt'
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

        if phot[i]=='IRAC2':
            ## DustPedia IRAC2 of M82 has no uncertainty map
            ## Create homogeneous uncertainty map
            bg = read_fits(path_tmp+src+'_'+phot[i]+'_DP').data[1720:1760,3110:3150]
            iuncert(path_tmp+src+'_'+phot[i]+'_DP',
                    filOUT=path_tmp+src+'_'+phot[i]+'_DP_unc', BG_image=bg)
        elif phot[i]=='IRAC3':
            Jy_per_pix_to_MJy_per_sr(raw_phot1+'_unc',
                                     filOUT=path_tmp+src+'_'+phot[i]+'_DP_unc')

        ## Reproject phot
        mtg = imontage('exact', tmpdir=path_tmp)
        mtg.reproject_mc(path_tmp+src+'_'+phot[i]+'_DP', refheader=coadd_fp,
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
        
        ## (Sect. 3.2) The weight maps contain the information on the number of frames 
        ## that were used to create the science mosaics at each pixel (value= # frames x10); 
        ## the pixel size of the weight maps is the same as the science mosaics.
        ## -- SINGS v5 release doc
        ## https://irsa.ipac.caltech.edu/data/SPITZER/SINGS/doc/sings_fifth_delivery_v2.pdf
        
        ## Only need to run ONCE for each photometry (comment after using)
        # gen_unc2 = input("Create SINGS IRAC2/IRAC3 uncertainty map (y/n)? ")
        # if gen_unc2=='y':
        #     ds_wgt = read_fits(wgt_phot2)
        #     if phot[i]== 'IRAC2':
        #         bg = read_fits(raw_phot2).data[542:582,1102:1142]
        #         bg_wgt = ds_wgt.data[542:582,1102:1142]
        #     elif phot[i]=='IRAC3':
        #         bg = read_fits(raw_phot2).data[402:442,1022:1062]
        #         bg_wgt = ds_wgt.data[402:442,1022:1062]
        #     wfac = .1
        #     iuncert(raw_phot2, filOUT=raw_phot2+'_unc',
        #             filWGT=wgt_phot2, wfac=wfac,
        #             BG_image=bg, BG_weight=bg_wgt)

        ## Reproject phot
        mtg = imontage('exact', tmpdir=path_tmp)
        mtg.reproject_mc(raw_phot2, refheader=coadd_fp,
                         dist='norm', Nmc=Nmc, filOUT=tmp_phot2)
    
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

    ## Fit correlation
    ##-----------------
    ## Data
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
        pix_phot2 = data_phot2[:,:].reshape((-1,))
        pix_spec = data_spec[:,:].reshape((-1,))
        pix_phot2_unc = unc_phot2[:,:].reshape((-1,))
        pix_spec_unc = unc_spec[:,:].reshape((-1,))
        
        # xgrid = np.arange(0,1e3,1)
        xgrid = np.logspace(-6,4,10000)
        
        ## NaNs mask
        mask = ~np.logical_or( np.isnan(pix_spec), np.isnan(pix_phot2) )
        
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
                   xlog=1, ylog=1, nonposx='clip', nonposy='clip',
                   xlim=(1e-3,1e4), ylim=(1e-3,1e4),
                   xlabel='DustPedia (MJy/sr)', ylabel='SINGS (MJy/sr)',
                   figsize=(8,8), title=src+'_'+phot[i],
                   titlesize=20, labelsize=10, ticksize=10, legendsize=10)
        p0.save(path_cal+'SINGS-DP_'+phot[i])
            
        ## MIR spec - SINGS
        p = pplot(pix_phot2, pix_spec,
                  yerr=pix_spec_unc, xerr=pix_phot2_unc,
                  fmt='.', c='k', ec='r', legend='upper left',
                  xlog=1, ylog=1, nonposx='clip', nonposy='clip',
                  xlim=(1e-6,1e4), ylim=(1e-6,1e4),
                  xlabel='SINGS (MJy/sr)', ylabel='MIR spec (MJy/sr)',
                  figsize=(8,8), title=src+'_'+phot[i],
                  titlesize=20, labelsize=10, ticksize=10, legendsize=10)

        ## Linear fit
        if i==0:
            # mask_irac2 = pix_spec/pix_phot2>1.e-2
            # mask = np.logical_and( mask, mask_irac3 )
            
            popt, pcov = curve_fit(f_lin0, pix_spec[mask], pix_phot2[mask],
                                   sigma=pix_phot2_unc[mask])
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
        else:
            # mask_irac3 = pix_spec/pix_phot2>1.e-2
            # mask = np.logical_and( mask, mask_irac3 )
            
            popt, pcov = curve_fit(f_lin0, pix_spec[mask], pix_phot2[mask],
                                   sigma=pix_phot2_unc[mask])
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
        
        p.save(path_cal+'intercalib_'+phot[i])

##----------------------------------------------------------

##                      Plot spectra

##----------------------------------------------------------
do_plot = input("Plot MIR spectra (y/n)? ")
if do_plot=='y':
    ds1 = read_fits(out_irc+'_0')
    ds2 = read_fits(out_irs+'_0')
    ds = read_fits(path_out+src+'_0', path_out+src+'_unc')
    Nw, Ny, Nx = ds.data.shape
    pp = intercalib(path_out+src+'_0')
    pp.read_filter(phot)
    wcen = pp.wcen
    for x in range(Nx):
        for y in range(Ny):
            if ~np.isnan(ds.data[:,y,x]).any():
            # if ~np.isnan(ds1.data[:,y,x]).any():
            #     if ~np.isnan(ds2.data[:,y,x]).any():
                p = pplot(ds.wave, ds.data[:,y,x], yerr=ds.unc[:,y,x],
                          c='k', lw=.7, ec='r', xlog=1, ylog=1,
                          legend='upper left', figsize=(12,10),
                          title=src+'_('+str(x)+','+str(y)+')',
                          titlesize=20, labelsize=10, ticksize=10, legendsize=10)
                p.add_plot(ds1.wave, ds1.data[:,y,x],
                           c='y', lw=.7, ls='-', zorder=-1, label='IRC raw')
                p.add_plot(ds2.wave, ds2.data[:,y,x],
                           c='y', lw=.7, ls='-', zorder=-1, label='IRS raw')
                
                p.add_plot(wcen[0], d1_phot2[y,x], yerr=u1_phot2[y,x],
                           c='m', marker='o', ms=10, zorder=100, label=phot[0])
                p.add_plot(wcen[0], d1_spec[y,x], yerr=u1_spec[y,x],
                           c='g', marker='^', ms=10, zorder=101, label='MIR spec-'+phot[0])
                p.add_plot(wcen[1], d2_phot2[y,x], yerr=u2_phot2[y,x],
                           c='orange', marker='o', ms=10, zorder=100, label=phot[1])
                p.add_plot(wcen[1], d2_spec[y,x], yerr=u2_spec[y,x],
                           c='c', marker='^', ms=10, zorder=101, label='MIR spec-'+phot[1])
                
                
                ph = ['IRAC1','IRAC4','MIPS1']
                pp.read_filter(ph)
                wc = pp.wcen
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
                    p.add_plot(wc[i], data_phot[y,x], yerr=unc_phot[y,x],
                               c=c1[i], marker='o', ms=10, zorder=100, label=ph[i])
                    if i==0:
                        p.add_plot(wc[i], data_spec[y,x], yerr=unc_spec[y,x],
                                   c=c2[i], marker='^', ms=10, zorder=101, label='IRC-'+ph[i])
                    else:
                        p.add_plot(wc[i], data_spec[y,x], yerr=unc_spec[y,x],
                                   c=c2[i], marker='^', ms=10, zorder=101, label='IRS-'+ph[i])
                
                xtic = [2, 3, 4, 5, 6, 7, 8, 10, 12, 15, 20, 30, 40]
                # xtic = []
                xtic_min = np.arange(2., 41., 1.)
                p.ax.set_xticks(xtic, minor=False) # major
                p.ax.set_xticks(xtic_min, minor=True) # minor
                p.ax.xaxis.set_major_formatter(ScalarFormatter()) # major
                p.ax.xaxis.set_minor_formatter(NullFormatter()) # minor
                
                p.save(path_fig+src+'_('+str(x+1)+','+str(y+1)+')', transparent=1)

##----------------------------------------------------------

##                     Final cut spectra

##----------------------------------------------------------
final_cut = input("MIR spectral map final cut with NaN edge, manuel mask, etc. (y/n)? ")
if final_cut=='y':
    if not os.path.exists(path_out+'final/'):
        os.makedirs(path_out+'final/')

    imp = improve(path_out+src+'_0')
    im = imp.im
    imp_unc = improve(path_out+src+'_unc')
    unc = imp_unc.im
    Nw, Ny, Nx = imp.im.shape
    for x in range(Nx):
        for y in range(Ny):
            if np.isnan(im[:,y,x]).any():
                imp.im[:,y,x] = np.nan
                imp_unc.im[:,y,x] = np.nan
    xlist = []
    for x in range(Nx):
        if not np.isnan(imp.im[:,:,x]).all():
            xlist.append(x)
    ylist = []
    for y in range(Ny):
        if not np.isnan(imp.im[:,y,:]).all():
            ylist.append(y)
    xmin = min(xlist)
    xmax = max(xlist)+1
    ymin = min(ylist)
    ymax = max(ylist)+1
    dx = xmax-xmin
    dy = ymax-ymin
    x0 = xmin+dx/2
    y0 = ymin+dy/2
    
    imp.crop(sizpix=(dx,dy), cenpix=(x0,y0))
    imp_unc.crop(sizpix=(dx,dy), cenpix=(x0,y0))

    ## Manual mask
    delcol = [1,5,9,11]
    for x in delcol:
        imp.im[:,:,x] = np.nan
        imp_unc.im[:,:,x] = np.nan
    imp.im[:3,:,:] = np.nan # AKARI wvl edge
    imp_unc.im[:3,:,:] = np.nan
    imp.im[225:265,:,:] = np.nan # AKARI-SL2 wvl edge
    imp_unc.im[225:265,:,:] = np.nan
    imp.im[:,:6,2] = np.nan
    imp_unc.im[:,:6,2] = np.nan
    imp.im[:,55:,2] = np.nan
    imp_unc.im[:,55:,2] = np.nan

    ## Print spectral map info
    Nw, Ny, Nx = imp.im.shape
    Npix = 0
    for y in range(Ny):
        for x in range(Nx):
            if not np.isnan(imp.im[:,y,x]).all():
                Npix += 1
    print('Total number of spectra: '+str(Npix))
    print('Size of spectral sampling: '+str(Nw))
    
    write_fits(path_out+'final/'+src, imp.hdr, imp.im, imp.wvl, imp.wmod)
    write_fits(path_out+'final/'+src+'_unc', imp.hdr, imp_unc.im, imp.wvl, imp.wmod)

    ds = read_fits(path_out+'final/'+src, path_out+'final/'+src+'_unc')

    Nw, Ny, Nx = ds.data.shape
    for x in range(Nx):
        for y in range(Ny):
            if ~np.isnan(ds.data[:,y,x]).all():
                p = pplot(ds.wave, ds.data[:,y,x], yerr=ds.unc[:,y,x],
                          c='k', lw=.7, ec='r', xlog=0, ylog=0,
                          xlim=(3.,22.), ylim=(-.1,np.nanmax(ds.data[:503,y,x])),
                          legend='upper left', figsize=(12,10),
                          title=src+'_('+str(x)+','+str(y)+')',
                          titlesize=20, labelsize=20, ticksize=20, legendsize=20)
                # p.add_plot(ds.wave, ds.data[:,y,x]/ds.unc[:,y,x], c='c')
                
                # xtic = [2, 3, 4, 5, 6, 7, 8, 10, 12, 15, 20, 30, 40]
                # xtic_min = np.arange(2., 41., 1.)
                xtic = np.arange(3., 21., 1)
                xtic_min = np.arange(3., 21., .1)
                p.ax.set_xticks(xtic, minor=False) # major
                p.ax.set_xticks(xtic_min, minor=True) # minor
                p.ax.xaxis.set_major_formatter(ScalarFormatter()) # major
                p.ax.xaxis.set_minor_formatter(NullFormatter()) # minor

                ## Line check
                cline = [
                         # 4.052,     # Bra
                         # 5.128657,  # HI6-10
                         # 5.51116,   # H2S7
                         # 5.6098,    # MgV1
                         # 5.908213,  # Huc
                         # 5.981,     # KIV
                         # 6.10856,   # H2S6
                         # 6.709,     # ClV
                         # 6.90952,   # H2S5
                         # 6.947984,  # HeII1
                         # 6.985274,  # ArII
                         # 7.3178,    # NaIII
                         # 7.459858,  # Pfa
                         # 7.502493,  # Hub
                         # 7.6524,    # NeVI
                         # 7.8145,    # FeVII
                         # 7.90158,   # ArV
                         # 8.02505,   # H2S4
                         # 8.760064,  # HI7-10
                         # 8.99103,   # ArIII1
                         # 9.042,     # NiVI
                         # 9.5267,    # FeVII
                         # 9.66491,   # H2S3
                         # 9.713475,  # HeII2
                         # 10.3385,   # SiI
                         # 10.5105,   # SIV
                         # 12.27861,  # H2S2
                         # 12.368527, # Hua
                         # 12.81355,  # NeII
                         # 13.10219,  # ArV
                         # 13.521,    # MgV2
                         # 14.32168,  # NeV1
                         # 15.555,    # NeIII1
                         # 17.0398,   # H2S1
                         # 17.608246, # HI11-18
                         # 18.7129,   # SIII1
                         # 19.061898, # HI7-8
                         # 21.8291,   # ArIII2
                         # 24.3175,   # NeV2
                         # 25.8903,   # OIV
                         # 25.9882,   # FeII1
                         # 28.21883,  # H2S0
                         # 33.4810,   # SIII2
                         # 34.8152,   # SiII
                         # 35.3491,   # FeII2
                         # 36.0135    # NeIII2
                         ]
                for clin in cline:
                    p.ax.axvline(clin, color='g')
                ## Band check
                cband = [
                         # 3.291,     # Main 3.3
                         # 3.399,     # Main 3.4
                         # 3.499,     # Small 3.5
                         # 5.2394667, # Small 5.2
                         5.6437461, # Small 5.7 (1)
                         5.7490305, # Small 5.7 (2)
                         6.0106598, # Small 6.0
                         6.2034273, # Main 6.2 (1)
                         6.2672596, # Main 6.2 (2)
                         6.6273788, # Small 6.6
                         6.8548833, # Small 6.8
                         7.0791725, # Small 7.1
                         7.6000000, # Plateau 7.7
                         7.6171123, # Main 7.7 (1)
                         7.8704769, # Main 7.7 (2)
                         8.3623706, # Small 8.3
                         8.6204540, # Main 8.6
                         # 9.5244838, # Small 9.5
                         10.707132, # Small 10.7
                         11.038349, # Small 11.0
                         11.237893, # Main 11.2
                         11.400432, # Plateau 11.3
                         # 11.796389, # Small 11.8
                         11.949674, # Small 11.9
                         12.626842, # Main 12.7 (1)
                         12.760273, # Main 12.7 (2)
                         13.559342, # Small 13.6
                         14.257133, # Small 14.2
                         # 15.893117, # Small 15.6
                         16.482868, # Small 16.4
                         17.082868, # Plateau 17.0
                         # 17.428485, # Small 17.4
                         17.771096, # Small 17.8
                         # 18.925630, # Small 18.9
                         ]
                for cban in cband:
                    p.ax.axvline(cban, color='b')
                
                p.save(path_out+'final/'+src+'_('+str(x+1)+','+str(y+1)+')', transparent=1)
