#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""

SOURCE: NGC 1569

"""

import os

## laputan
from laputan.astrom import fixwcs
from laputan.inout import fitsext
from laputan.imaging import iswarp


##----------------------------------------------------------

##                    Inputs & Outputs

##----------------------------------------------------------
src = 'NGC1569'
Nmc = 20
verbose = False

## Current dir
# croot = os.getcwd()+'/'
## Path of current file
croot = os.path.dirname(os.path.abspath(__file__))+'/'
## MIRAGE root path
mroot = '/Users/dhu/Github/MISSILE/MIRAGE/'

## IDL dir
##---------
path_idl = mroot+'idl/'

## Work dir
##----------
path_work = '/Users/dhu/Data/'

## Data dir
path_irc = path_work+'AKARI/'+src+'/'
path_irs = path_work+'Spitzer/'+src+'/'
path_phot = path_work+'Photometry/'+src+'/'
path_ker = path_work+'Kernels/' # idl/convolve_image.pro

## Tmp files
##-----------
path_tmp = path_work+'PAHPedia/tmp/'
if not os.path.exists(path_tmp):
    os.makedirs(path_tmp)
path_build = path_tmp+'cubuild/' # IRC cube building
if not os.path.exists(path_build):
    os.makedirs(path_build)
path_conv = path_tmp+'conv/' # idl/convolve_image.pro
if not os.path.exists(path_conv):
    os.makedirs(path_conv)
csv_ker = path_tmp+'kernelist' # idl/conv_prog.pro

## Outputs
##---------
path_out = path_work+'PAHPedia/'+src+'/'
if not os.path.exists(path_out):
    os.makedirs(path_out)

## Inter-calibrations
path_cal = path_out+'calib/'
if not os.path.exists(path_cal):
    os.makedirs(path_cal)

## Plots
path_fig = path_out+'figures/'
if not os.path.exists(path_fig):
    os.makedirs(path_fig)


##----------------------------------------------------------

##                         IRC data

##----------------------------------------------------------
out_irc = path_out+src+'_IRC'

slits = ['Ns', 'Nh']
obsid = [
    ['3390001.1','F011100297_N002','NG'], # A1-6, H1-2
    ['3390001.2','F011100338_N002','NG'], # 
    ['3390002.1','F011100795_N002','NG'], # C1-6, I1-2
    # ['3390002.2','F011176073_N002','NG'], # Matching failed
    ['3390003.1','F011100379_N002','NG'], # B1-6, G1-2
    ['5124077.1','F007174142_N002','NG'], # cap
    ['5125401.1','F010117172_N002','NG'], # J1-2, E1-3
    # ['5125402.1','F010117172_N002','NP'], # Matching failed
    ['5125403.1','F010116924_N002','NG'], # 
    # ['5125404.1','F010117338_N002','NP'], # Matching failed
    ['5125405.1','F010116950_N002','NG'], # D1-3, F1-2
    # ['5125406.1','F010117086_N002','NP'], # NP has 68 wvl instead of 259
]

fits_irc = []
parobs = []
for obs in obsid:
    for s in slits:
        fits_irc.append(path_build+obs[0]+'_'+s)
        parobs.append([ obs[0], s, obs[1], obs[2] ])


##----------------------------------------------------------

##                  IRS data (via CUBISM)

##----------------------------------------------------------
out_irs = path_out+src+'_IRS'

# chnl = ['SL2', 'SL1', 'LL2', 'LL1', 'SL3', 'LL3']
chnl = ['SH', 'LH']

lab_sl2 = []
lab_sl3 = []
lab_sl1 = []
lab_ll2 = []
lab_ll3 = []
lab_ll1 = []
lab_sh = []
lab_lh = []

fits_irs = []
for ch in chnl:
    if ch=='SL2':
        lab_ch = lab_sl2
    elif ch=='SL3':
        lab_ch = lab_sl3
    elif ch=='SL1':
        lab_ch = lab_sl1
    elif ch=='LL2':
        lab_ch = lab_ll2
    elif ch=='LL3':
        lab_ch = lab_ll3
    elif ch=='LL1':
        lab_ch = lab_ll1
    elif ch=='SH':
        lab_ch = lab_sh
    elif ch=='LH':
        lab_ch = lab_lh
        
    fits_ch = []
    if len(lab_ch)==0:
        fits_ch.append(path_irs+src+'_'+ch)
    else:
        for t in lab_ch:
            fits_ch.append(path_irs+src+'_'+t+'_'+ch)
            
    fits_irs.append(fits_ch)


##----------------------------------------------------------

##                   Data homogenisation

##----------------------------------------------------------
## Reprojection
##--------------
# coadd_tool = 'swarp'
coadd_tool = 'reproject'

## IRC grid
# coadd_fp = fixwcs(out_irc+fitsext).header
## IRS grid
refheader = fixwcs(fits_irs[0][0]+fitsext).header
swp = iswarp(sum(fits_irs, []), refheader=refheader,
             tmpdir=path_tmp, verbose=verbose)
coadd_fp = swp.refheader

## Convolution
##-------------
fits_ker = []
fwhm = [1., 1.5, 2., 2.5, 3., 3.5, 4., 4.5, 5., 5.5, 6.]
# psf_ref = 'IRAC_5.8' # 2.11 (< LL1)
# psf_ref = 'IRAC_8.0'# 2.82 (< LL1)
psf_ref = 'Gauss_06.0'
# psf_ref = 'MIPS_24' # 6.43
# psf_ref = 'WISE_ATLAS_11.6' # 6.60"
# psf_ref = 'WISE_ATLAS_22.1' # 11.89"
for w in fwhm:
    fits_ker.append(path_ker+'Kernel_HiRes_Gauss_0'+
                    str(w)+'_to_'+psf_ref)
