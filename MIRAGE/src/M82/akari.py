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
import math
import numpy as np
from scipy.optimize import curve_fit
## rapyuta
from rapyuta.inout import fclean, fitsext, read_fits, write_fits, read_hdf5
from rapyuta.imaging import ( iconvolve, iswarp, cupid, iuncert,
                              icrop, igroupixel, ismooth, 
                              imontage, improve, Jy_per_pix_to_MJy_per_sr )
from rapyuta.astrom import fixwcs
from rapyuta.arrays import closest, pix2sup, sup2pix
from rapyuta.calib import intercalib
from rapyuta.maths import f_lin, f_lin0, f_lin1
from rapyuta.plots import pplot


##----------------------------------------------------------

##                     Preliminaries

##----------------------------------------------------------

## Local
from buildinfo import ( src, Nmc, verbose, coadd_tool, colors, markers,
                        path_irc, path_build, parobs, # build
                        filog, fits_irc, out_irc, path_out,
                        path_phot, path_ker, path_cal, # calib
                        path_idl, csv_ker, path_tmp, path_fig
)

## Banner
print('\n============================================================\n')

print('        MIRAGE - AKARI/IRC cube builder - '+src)

print('\n============================================================\n')

Nobs = len(parobs)

# Nmc = 2

# exit()
    
##----------------------------------------------------------

##                      Build slits

##----------------------------------------------------------
do_build = input("(Re)build IRC slits (y/n)? ")
if do_build=='y':

    for i in trange(Nobs, #leave=False,
                    desc='<cupid> IRC slit building'):
        for j in trange(Nmc+1, leave=False,
                        desc=' - IRC slit unc [MC]'):
            cup = cupid(path_irc, obsid=parobs[i][1], slit=parobs[i][4],
                        spec=parobs[i][2], imref=parobs[i][3])
            if j==0:
                cup.spec_build(fits_irc[i]+'_'+str(j), tmpdir=path_build,
                               filRAW=fits_irc[i], fiLOG=filog+parobs[i][0],
                               Nx=parobs[i][7], Ny=parobs[i][5], Nsub=parobs[i][6],
                               pixscale=1, wmin=2.55, wmax=5.5, supix=True)
                if i==0:
                    wave0 = cup.wave()
            else:
                cup.spec_build(fits_irc[i]+'_'+str(j), tmpdir=path_build,
                               dist='splitnorm', sig_pt=3, fill_pt='med',
                               Nx=parobs[i][7], Ny=parobs[i][5], Nsub=parobs[i][6],
                               pixscale=1, wmin=2.55, wmax=5.5, supix=True)
            ## Interpolate wavelength grid of IRC slits
            ##------------------------------------------
            ## (A,B,C,D,G,H,I,J,L,N) (E,F,K) (M) have slightly different wgrid...
            ismooth(fits_irc[i]+'_'+str(j), wgrid=wave0, filOUT=fits_irc[i]+'_'+str(j))
            
        ## MC unc
        if Nmc>1:
            mcimage = []
            for j in range(Nmc):
                ds = read_fits(fits_irc[i]+'_'+str(j+1))
                mcimage.append(ds.data)
            mcimage = np.array(mcimage)
            unc = np.nanstd(mcimage, axis=0)
            write_fits(out_irc[i]+'_unc', ds.header, unc, ds.wave)


##----------------------------------------------------------

##                 Coadd slits (for footprint)

##----------------------------------------------------------
irc_coadd = input("Coadd IRC slits (y/n)? ")
if irc_coadd=='y':

    ## Coadding
    ##----------

    ## IRC grid
    refheader = fixwcs(fits_irc[0]+fitsext).header
    swp = iswarp(fits_irc, refheader=refheader, 
                 tmpdir=path_build, verbose=verbose)
    footprint_irc = swp.refheader

    slice_irc = []
    for f in fits_irc:
        ds = read_fits(f)
        ds.data[0] = 1
        write_fits(f+'_slice', ds.header, ds.data[0])
        slice_irc.append(f+'_slice')

    if coadd_tool=='swarp':
        
        ## <iswarp> coadding
        ##===================
        swp = iswarp(refheader=footprint_irc,
                     tmpdir=path_build, verbose=verbose)
        swp.combine(slice_irc,# combtype='wgt_avg',
                    keepedge=True, #cropedge=True,
                    filOUT=path_out+src+'_footprint')

    elif coadd_tool=='reproject':

        ## <imontage> coadding
        ##=====================
        mtg = imontage('exact', tmpdir=path_tmp, verbose=verbose)
        mtg.coadd(slice_irc, refheader=footprint_irc,
                  filOUT=path_out+src+'_footprint')

    ## Crop NaN edge
    edge = 2 # leave an edge of some pixels
    data = read_fits(path_out+src+'_footprint').data
    Ny, Nx = data.shape
    xlist = []
    for x in range(Nx):
        if not np.isnan(data[:,x]).all():
            xlist.append(x)
    ylist = []
    for y in range(Ny):
        if not np.isnan(data[y,:]).all():
            ylist.append(y)
    xmin = min(xlist) - edge
    if xmin<0:
        xmin = 0
    xmax = max(xlist)+1 + edge
    if xmax>Nx:
        xmax = Nx
    ymin = min(ylist) - edge
    if ymin<0:
        ymin = 0
    ymax = max(ylist)+1 + edge
    if ymax>Ny:
        ymax = Ny
    dx = xmax-xmin
    dy = ymax-ymin
    x0 = xmin+dx/2
    y0 = ymin+dy/2
    
    icrop(path_out+src+'_footprint', filOUT=path_out+src+'_footprint',
          sizpix=(dx,dy), cenpix=(x0,y0), verbose=verbose)
    
## Header of th atlas IRC map
header_atlas = fixwcs(path_out+src+'_footprint'+fitsext).header


##----------------------------------------------------------

##                   Inter-calibration

##----------------------------------------------------------
pre_calib1 = input("Prepare photometry: DustPedia IRAC1 (y/n)? ")
pre_calib2 = input("Prepare photometry: SINGS IRAC1 (y/n)? ")
synt_phot = input("Prepare IRC synthetic photometry (y/n)? ")
do_intcal = input("Do inter-calibration (y/n)? ")
if (synt_phot!='y' or pre_calib2!='y') and do_intcal=='y':
    warnings.warn('You need to prepare (synthetic) photometry if you have not done it yet.')

os.makedirs(path_tmp+'calib/', exist_ok=True)
    
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
    bg = read_fits(tmp_phot1).data[1713:1753,3110:3150]*10
    iuncert(tmp_phot1, filOUT=tmp_phot1+'_unc', BG_image=bg)

    ## Reproject phot (atlas)
    if coadd_tool=='swarp':
        swp = iswarp(refheader=header_atlas, tmpdir=path_tmp)
        swp.combine_mc(tmp_phot1, keepedge=True, #cropedge=True,
                       dist='norm', Nmc=Nmc, filOUT=tmp_phot1+'_MC')
    elif coadd_tool=='reproject':
        mtg = imontage('exact', tmpdir=path_tmp)
        mtg.reproject_mc(tmp_phot1, refheader=header_atlas,
                         dist='norm', Nmc=Nmc, filOUT=tmp_phot1+'_MC')
    
    for i in trange(Nobs, #leave=False,
                    desc=phot+' (DustPedia) processing'):
        xscale, yscale = read_hdf5(filog+parobs[i][0], 'Super pixel size')
        refheader = fixwcs(fits_irc[i]+'_0'+fitsext).header

        for j in trange(Nmc+1, leave=False,
                        desc=' - '+parobs[i][0]+' [MC]'):
            if i==0:
                ## Convolve phot
                if j==0:
                    conv = iconvolve(tmp_phot1+'_MC',
                                     kfile=phot_ker, klist=csv_ker,
                                     filOUT=tmp_phot1+'_MC')
                else:
                    conv = iconvolve(tmp_phot1+'_MC_'+str(j),
                                     kfile=phot_ker, klist=csv_ker,
                                     # filUNC=tmp_phot1+'_unc', dist='norm',
                                     filOUT=tmp_phot1+'_MC_'+str(j))
                conv.do_conv(idldir=path_idl)

            ## Reproject phot (slit)
            if coadd_tool=='swarp':
                swp = iswarp(refheader=refheader, tmpdir=path_tmp)
                if j==0:
                    swp.combine(tmp_phot1+'_MC',
                                filOUT=fits_phot1+'_'+parobs[i][0])
                else:
                    swp.combine(tmp_phot1+'_MC_'+str(j),
                                filOUT=tmp_phot1+'_'+parobs[i][0]+'_'+str(j))
            elif coadd_tool=='reproject':
                mtg = imontage('exact', tmpdir=path_tmp)
                if j==0:
                    mtg.reproject(tmp_phot1+'_MC', refheader=refheader,
                                  filOUT=fits_phot1+'_'+parobs[i][0])
                else:
                    mtg.reproject(tmp_phot1+'_MC_'+str(j), refheader=refheader,
                                  filOUT=tmp_phot1+'_'+parobs[i][0]+'_'+str(j))
            
            ## Calculate unc
            if j==0:
                igroupixel(fits_phot1+'_'+parobs[i][0],
                           xscale=xscale, yscale=yscale,
                           filOUT=fits_phot1+'_'+parobs[i][0])
            else:
                igroupixel(tmp_phot1+'_'+parobs[i][0]+'_'+str(j),
                           xscale=xscale, yscale=yscale,
                           filOUT=tmp_phot1+'_'+parobs[i][0]+'_'+str(j))
        if Nmc>1:
            mcimage = []
            for j in range(Nmc):
                ds = read_fits(tmp_phot1+'_'+parobs[i][0]+'_'+str(j+1))
                mcimage.append(ds.data)
            mcimage = np.array(mcimage)
            unc = np.nanstd(mcimage, axis=0)
            write_fits(fits_phot1+'_'+parobs[i][0]+'_unc', ds.header, unc)

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

    ## Reproject phot (atlas)
    if coadd_tool=='swarp':
        swp = iswarp(refheader=header_atlas, tmpdir=path_tmp)
        swp.combine_mc(raw_phot2, keepedge=True, #cropedge=True,
                       dist='norm', Nmc=Nmc, filOUT=tmp_phot2+'_MC')
    elif coadd_tool=='reproject':
        mtg = imontage('exact', tmpdir=path_tmp)
        mtg.reproject_mc(raw_phot2, refheader=header_atlas,
                         dist='norm', Nmc=Nmc, filOUT=tmp_phot2+'_MC')

    for i in trange(Nobs, #leave=False,
                    desc=phot+' (SINGS) processing'):
        xscale, yscale = read_hdf5(filog+parobs[i][0], 'Super pixel size')
        refheader = fixwcs(fits_irc[i]+'_0'+fitsext).header

        for j in trange(Nmc+1, leave=False,
                        desc=' - '+parobs[i][0]+' [MC]'):
            if i==0:
                ## Convolve phot
                if j==0:
                    conv = iconvolve(tmp_phot2+'_MC',
                                     kfile=phot_ker, klist=csv_ker,
                                     filOUT=tmp_phot2+'_MC')
                else:
                    conv = iconvolve(tmp_phot2+'_MC_'+str(j),
                                     kfile=phot_ker, klist=csv_ker,
                                     # filUNC=tmp_phot2+'_unc', dist='norm',
                                     filOUT=tmp_phot2+'_MC_'+str(j))
                conv.do_conv(idldir=path_idl)

            ## Reproject phot (slit)
            if coadd_tool=='swarp':
                swp = iswarp(refheader=refheader, tmpdir=path_tmp)
                if j==0:
                    swp.combine(tmp_phot2+'_MC',
                                filOUT=fits_phot2+'_'+parobs[i][0])
                else:
                    swp.combine(tmp_phot2+'_MC_'+str(j),
                                filOUT=tmp_phot2+'_'+parobs[i][0]+'_'+str(j))
            elif coadd_tool=='reproject':
                mtg = imontage('exact', tmpdir=path_tmp)
                if j==0:
                    mtg.reproject(tmp_phot2+'_MC', refheader=refheader,
                                  filOUT=fits_phot2+'_'+parobs[i][0])
                else:
                    mtg.reproject(tmp_phot2+'_MC_'+str(j), refheader=refheader,
                                  filOUT=tmp_phot2+'_'+parobs[i][0]+'_'+str(j))
            
            ## Calculate unc
            if j==0:
                igroupixel(fits_phot2+'_'+parobs[i][0],
                           xscale=xscale, yscale=yscale,
                           filOUT=fits_phot2+'_'+parobs[i][0])
            else:
                igroupixel(tmp_phot2+'_'+parobs[i][0]+'_'+str(j),
                           xscale=xscale, yscale=yscale,
                           filOUT=tmp_phot2+'_'+parobs[i][0]+'_'+str(j))
        if Nmc>1:
            mcimage = []
            for j in range(Nmc):
                ds = read_fits(tmp_phot2+'_'+parobs[i][0]+'_'+str(j+1))
                mcimage.append(ds.data)
            mcimage = np.array(mcimage)
            unc = np.nanstd(mcimage, axis=0)
            write_fits(fits_phot2+'_'+parobs[i][0]+'_unc', ds.header, unc)

## Synthetic photometry (spec)
##=============================
if (synt_phot=='y'):# and do_build=='y'):
    ## fits_irc[i] can be changed by inter-calib,
    ## One need to rebuild before synthetic photometry
    for i in trange(Nobs, #leave=False,
                    desc=phot+' (IRC synt phot) slits'):
        xscale, yscale = read_hdf5(filog+parobs[i][0], 'Super pixel size')
        for j in trange(Nmc+1, leave=False,
                        desc=' - '+parobs[i][0]+': IRC sythetic photometry [MC]'):
            ic = intercalib(fits_irc[i]+'_'+str(j))
            sp = ic.synthetic_photometry(phot, xscale=xscale, yscale=yscale)
            if j==0:
                write_fits(fits_spec+'_'+parobs[i][0], ic.hdr, sp.Fnu_filt)
            else:
                write_fits(tmp_spec+'_'+parobs[i][0]+'_'+str(j), ic.hdr, sp.Fnu_filt)
        if Nmc>1:
            mcimage = []
            for j in range(Nmc):
                ds = read_fits(tmp_spec+'_'+parobs[i][0]+'_'+str(j+1))
                mcimage.append(ds.data)
            mcimage = np.array(mcimage)
            unc = np.nanstd(mcimage, axis=0)
            write_fits(fits_spec+'_'+parobs[i][0]+'_unc', ds.header, unc)

##-----------------
## Fit correlation
##-----------------
## Data dictionary
## e.g. [{'A0': a0, 'A1': a1, ...}, {'B0': b0, 'B1': b1, ...}, ...]
dict_phot1 = []
dict_phot1_unc = []
dict_phot2 = []
dict_phot2_unc = []
dict_spec = []
dict_spec_unc = []
## Data 1D array
## e.g. [array([a1, a2, ...]), array([b1, b2, ...]), ...]
pix_phot1 = []
pix_phot1_unc = []
pix_phot2 = []
pix_phot2_unc = []
pix_spec = []
pix_spec_unc = []

# xgrid = np.arange(0,1e3,1)
xgrid = np.logspace(-2,3,1000)

for i in range(Nobs):
    # Nx = read_hdf5(filog+parobs[i][0], 'Slit width')
    Ny = read_hdf5(filog+parobs[i][0], 'Slit length')[0]
    xscale, yscale = read_hdf5(filog+parobs[i][0], 'Super pixel size')
    # Nxs = math.ceil(Nx/xscale)
    Nys = math.ceil(Ny/yscale)

    ## Init dict
    data1 = {}
    unc1 = {}
    data2 = {}
    unc2 = {}
    data3 = {}
    unc3 = {}
    for ys in range(Nys):
        subname = parobs[i][0][0]+str(ys+1)
        yarr = sup2pix(ys, yscale, Npix=Ny, origin=0)
        ## Phot1
        ds = read_fits(fits_phot1+'_'+parobs[i][0], fits_phot1+'_'+parobs[i][0]+'_unc')
        data1[subname] = ds.data[yarr[0],0]
        unc1[subname] = ds.unc[yarr[0],0]
        ## Phot2
        ds = read_fits(fits_phot2+'_'+parobs[i][0], fits_phot2+'_'+parobs[i][0]+'_unc')
        data2[subname] = ds.data[yarr[0],0]
        unc2[subname] = ds.unc[yarr[0],0]
        ## spec
        ds = read_fits(fits_spec+'_'+parobs[i][0], fits_spec+'_'+parobs[i][0]+'_unc')
        data3[subname] = ds.data[yarr[0],0]
        unc3[subname] = ds.unc[yarr[0],0]
    dict_phot1.append(data1)
    dict_phot1_unc.append(unc1)
    dict_phot2.append(data2)
    dict_phot2_unc.append(unc2)
    dict_spec.append(data3)
    dict_spec_unc.append(unc3)

    ## Dictionary to 1D array
    pix_phot1.append(np.array(list(dict_phot1[i].values())))
    pix_phot1_unc.append(np.array(list(dict_phot1_unc[i].values())))
    pix_phot2.append(np.array(list(dict_phot2[i].values())))
    pix_phot2_unc.append(np.array(list(dict_phot2_unc[i].values())))
    pix_spec.append(np.array(list(dict_spec[i].values())))
    pix_spec_unc.append(np.array(list(dict_spec_unc[i].values())))

    ## Mask NaNs and negtive values
    mask_nan = ~np.logical_or( np.isnan(pix_spec[i]), np.isnan(pix_phot2[i]) )
    mask_neg = np.logical_and( pix_spec[i]>0, pix_phot2[i]>0 )
    mask = np.logical_and( mask_nan, mask_neg )

    ## Add calibration error
    pix_phot1_unc[i] = np.sqrt( pix_phot1_unc[i]**2 + (pix_phot1[i]*.03)**2 ) # IRAC 3% (Carey2010)
    pix_phot2_unc[i] = np.sqrt( pix_phot2_unc[i]**2 + (pix_phot2[i]*.03)**2 )
    pix_spec_unc[i] = np.sqrt( pix_spec_unc[i]**2 + (pix_spec[i]*.2)**2 ) # IRC 10-20% (Ohyama+2007)
    
    ## S/N ratio
    # print(parobs[i][0]+' S/N (IRC) = \n', pix_spec[i][mask]/pix_spec_unc[i][mask])
    # print(parobs[i][0]+' S/N (SINGS) = \n', pix_phot2[i][mask]/pix_phot2_unc[i][mask])

    if mask.any():
        ## DP - SINGS plot
        if i==0:
            p0 = pplot(fmt='s', ec='grey', elw=1,# clib=colors,
                       xlog=1, ylog=1, nonposx='clip', nonposy='clip', tkform='mylog',
                       # xlim=(1e-2,1e3), ylim=(1e-2,1e3),
                       xlabel='DustPedia (MJy/sr)', ylabel='SINGS (MJy/sr)',
                       # title=src+' '+phot+' calibration',
                       figsize=(8,8),# right=.8, left=.15, bottom=.15,
                       loc='upper left', anchor=(0,1), legendalpha=0,
                       titlesize=20, xysize=20, tksize=20, legendsize=15)
        p0.add_plot(pix_phot1[i][mask], pix_phot2[i][mask],
                    yerr=pix_phot2_unc[i][mask], xerr=pix_phot1_unc[i][mask],
                    fmt='s', ec='grey', c=colors[i+1],
                    marker=markers[i], markersize=20, capsize=2,
                    label=parobs[i][0])
        p0.add_plot(pix_phot1[i][mask], pix_phot2[i][mask],
                    yerr=pix_phot2_unc[i][mask], xerr=pix_phot1_unc[i][mask],
                    fmt='o', ec='grey', c=colors[i+1], zorder=100,
                    marker=markers[i], markersize=.1, capsize=2) # put errorbar on the front layer
        p0.save(path_cal+'DP-SINGS_'+phot, transparent=True, figtight=True)
    
        if do_intcal=='y':
            
            ## IRC - SINGS plot
            if i==0:
                p = pplot(fmt='s', ec='grey', elw=1,# clib=colors,
                          xlog=1, ylog=1, nonposx='clip', nonposy='clip', tkform='mylog',
                          # xlim=(1e-2,1e3), ylim=(1e-2,1e3),
                          # title=src+' IRC-'+phot+' inter-calibration',
                          xlabel=r'$\rm IRC-sIRAC_{3.6\mu m}\ (MJy/sr)$',
                          ylabel=r'$\rm IRAC_{3.6\mu m}\ (MJy/sr)$',
                          figsize=(8,8),# right=.78, left=.12, bottom=.1, top=.95,
                          loc='upper left', anchor=(0,1), legendalpha=0,
                          titlesize=20, xysize=20, tksize=20, legendsize=15)
                                
            p.add_plot(pix_spec[i][mask], pix_phot2[i][mask],
                       yerr=pix_phot2_unc[i][mask], xerr=pix_spec_unc[i][mask],
                       fmt='s', ec='grey', c=colors[i+1],
                       marker=markers[i], markersize=20, capsize=2,
                       label=parobs[i][0])
            p.add_plot(pix_spec[i][mask], pix_phot2[i][mask],
                       yerr=pix_phot2_unc[i][mask], xerr=pix_spec_unc[i][mask],
                       fmt='s', ec='grey', c=colors[i+1], zorder=100,
                       marker=markers[i], markersize=.1, capsize=2) # put errorbar on the front layer
            p.save(path_cal+'IC_'+phot+'.png', transparent=True, figtight=True)
    
        ## Linear fit (IRC - SINGS, SLITS)
        ##=================================
        popt, pcov = curve_fit(f_lin0, pix_spec[i][mask], pix_phot2[i][mask],)
                               # sigma=pix_phot2_unc[i][mask])
        slit_gain = popt[0]
        # slit_gain = 1.
        # slit_off = popt[1]
        slit_off = 0.
        # slit_off = popt[0]
        
        if do_intcal=='y':
            print(parobs[i][0]+' inter-calibration ('+phot+') gain = {:.4}'.format(slit_gain))
            # print(parobs[i][0]+' inter-Calibration ('+phot+') offset = {:.4}'.format(slit_off))
            # label = parobs[i][0]+': y={0:.4}x'.format(slit_gain)
            # label = parobs[i][0]+': y={0:.4}x+{1:.4}'.format(slit_gain, slit_off)
            # p.add_plot(xgrid, f_lin0(xgrid, *popt),
            #            c='k', ls='-', label=label)
            # p.save(path_cal+'IC_'+phot+'.png', transparent=True, figtight=True)

## Linear fit (DP - SINGS)
##=========================
fitx = np.concatenate(pix_phot1, axis=0)
uncx = np.concatenate(pix_phot1_unc, axis=0)
fity = np.concatenate(pix_phot2, axis=0)
uncy = np.concatenate(pix_phot2_unc, axis=0)
maskfit = ~np.logical_or( np.isnan(fitx), np.isnan(fity) )

popt, pcov = curve_fit(f_lin0, fitx[maskfit], fity[maskfit],)
                       # sigma=uncy[maskfit])
atlas_gain = popt[0]
# atlas_gain = 1.
# atlas_off = popt[1]
atlas_off = 0.
# atlas_off = popt[0]

print('DustPedia - SINGS calibration ('+phot+') gain = {:.4}'.format(atlas_gain))
# print('DustPedia - SINGS calibration ('+phot+') offset = {:.4}'.format(atlas_off))
label = 'y={0:.4}x'.format(atlas_gain)
# label = 'y={0:.4}x+{1:.4}'.format(atlas_gain, atlas_off)
p0.add_plot(xgrid, f_lin0(xgrid, *popt),
            c='k', ls='-', label=label)
p0.ax.legend(loc='upper left', bbox_to_anchor=(1,1),
             fontsize=20, framealpha=0,)
p0.save(path_cal+'DP-SINGS_'+phot+'.png')

## Linear fit (IRC - SINGS, ATLAS)
##=================================
fitx = np.concatenate(pix_spec, axis=0)
uncx = np.concatenate(pix_spec_unc, axis=0)
fity = np.concatenate(pix_phot2, axis=0)
uncy = np.concatenate(pix_phot2_unc, axis=0)
maskfit = ~np.logical_or( np.isnan(fitx), np.isnan(fity) )

popt, pcov = curve_fit(f_lin0, fitx[maskfit], fity[maskfit],)
                       # sigma=uncy[maskfit])
atlas_gain = popt[0]
# atlas_gain = 1.
# atlas_off = popt[1]
atlas_off = 0.
# atlas_off = popt[0]

if do_intcal=='y':
    print('Atlas inter-calibration ('+phot+') gain = {:.4}'.format(atlas_gain))
    # print('Atlas inter-Calibration ('+phot+') offset = {:.4}'.format(atlas_off))
    label = 'y={0:.4}x'.format(atlas_gain)
    # label = 'y={0:.4}x+{1:.4}'.format(atlas_gain, atlas_off)
    p.add_plot(xgrid, f_lin0(xgrid, *popt),
               c='k', ls='-', label=label)
    p.ax.text(.9,.05,'(b)',size=30,c='grey',transform=p.ax.transAxes) # for the use of Hu_thesis
    p.ax.legend(loc='upper left', bbox_to_anchor=(1,1),
                fontsize=20, framealpha=0,)
    p.save(path_cal+'IC_'+phot+'.png', transparent=True, figtight=True)

## Spectral correction
##---------------------
write_irc = input("Write IRC spectra (y/n)? ")
if write_irc=='y':
    if do_build=='y':
        correct_atlas = input(" - ATLAS level correction (y/n)? ")
        if correct_atlas=='y':
            correct_pixel = None
        else:
            correct_pixel = input(" - PIXEL level correction (y/n)? ")
    else:
        warnings.warn('Rebuild slits before spectral correction.')
        correct_atlas = None
        correct_pixel = None
else:
    correct_atlas = None
    correct_pixel = None

for i in range(Nobs):
    xscale, yscale = read_hdf5(filog+parobs[i][0], 'Super pixel size')
    mask_nan = ~np.logical_or( np.isnan(pix_spec[i]), np.isnan(pix_phot2[i]) )
    mask_neg = np.logical_and( pix_spec[i]>0, pix_phot2[i]>0 )
    mask = np.logical_and( mask_nan, mask_neg )
    sc = intercalib(fits_irc[i]+'_0') # Attention: NOT fits_irc[i]
    Nw, Ny, Nx = sc.im.shape
    pixel_gain = np.ones((Ny,Nx))
    if mask.any():
        for y in range(Ny):
            ys = pix2sup(y, yscale, origin=0)
            if mask[ys]:
                pixel_gain[y,:] *= pix_phot2[i][ys]/pix_spec[i][ys]
            # if y%yscale==0:
            #     print(parobs[i][0]+' inter-calibration gain: {}'.format(pixel_gain[y,0]))
    if write_irc=='y':
        if correct_pixel=='y':
            calib_gain = pixel_gain
        elif correct_atlas=='y':
            calib_gain = atlas_gain
        else:
            calib_gain = 1.0
            print('No correction!')
        sc.correct_spec(calib_gain, filOUT=out_irc[i])


##----------------------
## Fit correlation (MC)
##----------------------
if do_intcal!='y':
    print('No correction!')
    
for j in trange(Nmc, #leave=False,
                desc='IRC spectral correction [MC]'):
    ## Data dictionary
    ## e.g. [{'A0': a0, 'A1': a1, ...}, {'B0': b0, 'B1': b1, ...}, ...]
    dict_spec = []
    dict_spec_unc = []
    ## Data 1D array
    ## e.g. [array([a1, a2, ...]), array([b1, b2, ...]), ...]
    pix_spec = []
    pix_spec_unc = []
    
    for i in range(Nobs):
        # Nx = read_hdf5(filog+parobs[i][0], 'Slit width')
        Ny = read_hdf5(filog+parobs[i][0], 'Slit length')[0]
        xscale, yscale = read_hdf5(filog+parobs[i][0], 'Super pixel size')
        # Nxs = math.ceil(Nx/xscale)
        Nys = math.ceil(Ny/yscale)
    
        ## Init dict
        data3 = {}
        unc3 = {}
        for ys in range(Nys):
            subname = parobs[i][0][0]+str(ys+1)
            yarr = sup2pix(ys, yscale, Npix=Ny, origin=0)
            ## spec
            ds = read_fits(tmp_spec+'_'+parobs[i][0]+'_'+str(j+1),
                           fits_spec+'_'+parobs[i][0]+'_unc') # intermediate unc
            data3[subname] = ds.data[yarr[0],0]
            unc3[subname] = ds.unc[yarr[0],0]
        dict_spec.append(data3)
        dict_spec_unc.append(unc3)
    
        ## Dictionary to 1D array
        pix_spec.append(np.array(list(dict_spec[i].values())))
        pix_spec_unc.append(np.array(list(dict_spec_unc[i].values())))
    
        ## Mask NaNs and negtive values
        mask_nan = ~np.logical_or( np.isnan(pix_spec[i]), np.isnan(pix_phot2[i]) )
        mask_neg = np.logical_and( pix_spec[i]>0, pix_phot2[i]>0 )
        mask = np.logical_and( mask_nan, mask_neg )
        
    ## Linear fit (IRC - SINGS, ATLAS)
    ##=================================
    fitx = np.concatenate(pix_spec, axis=0)
    uncx = np.concatenate(pix_spec_unc, axis=0)
    fity = np.concatenate(pix_phot2, axis=0)
    uncy = np.concatenate(pix_phot2_unc, axis=0)
    maskfit = ~np.logical_or( np.isnan(fitx), np.isnan(fity) )
    
    popt, pcov = curve_fit(f_lin0, fitx[maskfit], fity[maskfit],)
                           # sigma=uncy[maskfit])
    atlas_gain = popt[0]
    # atlas_gain = 1.
    # atlas_off = popt[1]
    atlas_off = 0.
    # atlas_off = popt[0]
    
    ## Spectral correction
    ##---------------------
    for i in range(Nobs):
        xscale, yscale = read_hdf5(filog+parobs[i][0], 'Super pixel size')
        mask_nan = ~np.logical_or( np.isnan(pix_spec[i]), np.isnan(pix_phot2[i]) )
        mask_neg = np.logical_and( pix_spec[i]>0, pix_phot2[i]>0 )
        mask = np.logical_and( mask_nan, mask_neg )
        sc = intercalib(fits_irc[i]+'_'+str(j+1))
        Nw, Ny, Nx = sc.im.shape
        pixel_gain = np.ones((Ny,Nx))
        if mask.any():
            for y in range(Ny):
                ys = pix2sup(y, yscale, origin=0)
                if mask[ys]:
                    pixel_gain[y,:] *= pix_phot2[i][ys]/pix_spec[i][ys]
                # if y%yscale==0:
                #     print(parobs[i][0]+' inter-calibration gain: {}'.format(pixel_gain[y,0]))
        if (write_irc=='y'):
            ## fits_irc[i] can be changed by inter-calib,
            ## One need to rebuild before synthetic photometry
            if correct_pixel=='y':
                calib_gain = pixel_gain
            elif correct_atlas=='y':
                calib_gain = atlas_gain
            else:
                calib_gain = 1.0
            sc.correct_spec(calib_gain, filOUT=fits_irc[i]+'_'+str(j+1))

if (write_irc=='y' and Nmc>1):
    for i in trange(Nobs, leave=False,
                    desc='Calculating uncertainties for IRC slits'):
        mcimage = []
        for j in range(Nmc):
            ds = read_fits(fits_irc[i]+'_'+str(j+1))
            mcimage.append(ds.data)
        mcimage = np.array(mcimage)
        unc = np.nanstd(mcimage, axis=0)
        write_fits(out_irc[i]+'_unc', ds.header, unc, ds.wave)


##----------------------------------------------------------

##                      Plot spectra

##----------------------------------------------------------
plot_spec = input("Plot IRC spectra (y/n)? ")
if plot_spec=='y':

    for i in range(Nobs):
        ds = read_fits(out_irc[i], out_irc[i]+'_unc')
        Nw, Ny, Nx = ds.data.shape
        ds0 = read_fits(fits_irc[i]+'_0')
        pp = intercalib(out_irc[i])
        pp.read_filter(phot)
        wcen = pp.wcen
        
        xscale, yscale = read_hdf5(filog+parobs[i][0], 'Super pixel size')
        for y in range(Ny):
            ys = pix2sup(y, yscale, origin=0)
            subname = parobs[i][0][0]+str(ys+1)
            maskspec = 1
            # maskspec = ~np.isnan(ds.data[:,y,0]).any()
            if (maskspec and y%yscale==0):
                xtic = [2.5, 3, 3.5, 4, 4.5, 5]
                xtic_min = np.arange(2.4, 5.2, .1)

                ## Add calibration error
                ds.unc[:,y,0] = np.sqrt( ds.unc[:,y,0]**2 + (ds.data[:,y,0]*.2)**2 ) # IRC 10-20% (Ohyama+2007)

                ## Individual spectra
                ##--------------------
                p0 = pplot(ds.wave, ds.data[:,y,0], yerr=ds.unc[:,y,0],
                           # xlog=1, ylog=1, tkform='mylog',
                           c='k', lw=2, ec='r', elw=1, label=subname,
                           xlabel=r'$\rm Wavelengths\ \lambda\ (\mu m)$',
                           ylabel=r'$\rm F_{\nu}\ (MJy/sr)$',
                           figsize=(12,8), loc='upper left', legendalpha=0,
                           # title='IRC_'+subname,
                           titlesize=20, xysize=20, tksize=20, legendsize=20)
                ## Non-intercalib
                p0.add_plot(ds.wave, ds0.data[:,y,0],
                            c='y', lw=2, ls='--', zorder=99)
                ## Photometry
                p0.add_plot(wcen[0], dict_phot2[i][subname], yerr=dict_phot2_unc[i][subname],
                            c='m', marker='o', ms=10, zorder=100, label=phot)
                ## Synthetic photometry
                p0.add_plot(wcen[0], dict_spec[i][subname], yerr=dict_spec_unc[i][subname],
                            c='g', marker='^', ms=10, zorder=101, label='IRC-'+phot)
                
                p0.save(path_fig+'IRC_'+subname, transparent=True, figtight=True)

                ## Overlapped in one figure
                ##--------------------------
                legendloc = 'upper left'
                if parobs[i][0]=='E1_4':
                    legendloc = 'upper right'
                if y==0:
                    p = pplot(ds.wave, ds.data[:,y,0], yerr=ds.unc[:,y,0],
                              # xlog=1, ylog=1, tkform='mylog',
                              lw=2, ec='grey', elw=1, label=subname,
                              clib=['c','m','y','r','orange','b','g'],
                              xlabel=r'$\rm Wavelengths\ \lambda\ (\mu m)$',
                              ylabel=r'$\rm F_{\nu}\ (MJy/sr)$',
                              xtk=xtic, xtkmi=xtic_min,# tkform='mylog',
                              figsize=(14,8), loc=legendloc, legendalpha=0,
                              titlesize=20, xysize=20, tksize=20, legendsize=20)
                else:
                    p.add_plot(ds.wave, ds.data[:,y,0], yerr=ds.unc[:,y,0],
                               lw=2, ec='grey', elw=1, label=subname)

                # p.ax.legend(loc='upper left',# bbox_to_anchor=(0,1),
                #             fontsize=20, framealpha=0)
                
                p.save(path_fig+src+'_IRC_'+subname[0], transparent=True, figtight=True)
