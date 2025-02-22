#!/usr/bin/env python
# -*- coding: utf-8 -*-

### DON NOT FORGET TO MODIFY idlib

import os

verbose = False

src = 'M82'
Nmc = 100

## Current dir
##-------------
path_cur = os.getcwd()+'/'
path_par = os.path.dirname(os.path.abspath(__file__))+'/' # param file path

## IDL dir
##---------
path_idl = path_par+'idlib/'

## Root dir
##----------
## Root of data and outputs
path_root = '/Users/dhu/Data/'
## Data dir
path_irs = path_root+'Spitzer/'+src+'/'
path_phot = path_root+'Photometry/'+src+'/'
path_ker = path_root+'Kernels/'

## IRS data (via CUBISM)
##-----------------------
fits_sl2 = []
fits_sl3 = []
fits_sl1 = []
fits_ll2 = []
fits_ll3 = []
fits_ll1 = []
fits_irs = []

chnl = ['SL2','SL1','LL2','LL1']
# chnl = ['SL2', 'SL3', 'SL1', 'LL2', 'LL3', 'LL1']
Nch = len(chnl) # Number of chnl (ordered as above) used
lab_sl = ['04','06S','06N','08','08c','09N3','09N2']
lab_ll = ['03','04','05','06','08','09N3','09N5','09N2']
# lab_sl = ['']
# lab_ll = ['']
for ch in chnl:
    ## SL
    if ch[:2]=='SL':
        for t in lab_sl:
            f = path_irs+src+'_'+t+'_'+ch
            if ch=='SL2':
                fits_sl2.append(f)
            if ch=='SL3':
                fits_sl3.append(f)
            if ch=='SL1':
                fits_sl1.append(f)
    ## LL
    elif ch[:2]=='LL':
        for t in lab_ll:
            f = path_irs+src+'_'+t+'_'+ch
            if ch=='LL2':
                fits_ll2.append(f)
            if ch=='LL3':
                fits_ll3.append(f)
            if ch=='LL1':
                fits_ll1.append(f)

if 'SL2' in chnl:
    fits_irs.append(fits_sl2)
if 'SL3' in chnl:
    fits_irs.append(fits_sl3)
if 'SL1' in chnl:
    fits_irs.append(fits_sl1)
if 'LL2' in chnl:
    fits_irs.append(fits_ll2)
if 'LL3' in chnl:
    fits_irs.append(fits_ll3)
if 'LL1' in chnl:
    fits_irs.append(fits_ll1)

## Convolution
##-------------
fits_ker = []
psf = [1., 1.5, 2., 2.5, 3., 3.5, 4., 4.5, 5., 5.5, 6.]
# psf_ref = 'IRAC_5.8' # 2.11 (< LL1)
# psf_ref = 'IRAC_8.0'# 2.82 (< LL1)
psf_ref = 'Gauss_06.0'
# psf_ref = 'MIPS_24' # 6.43"
# psf_ref = 'WISE_MAP_11.6' # 6.60"
# psf_ref = 'WISE_MAP_22.1' # 11.89"
for p in psf:
    fits_ker.append(path_ker+'Kernel_HiRes_Gauss_0'+
                    str(p)+'_to_'+psf_ref)

## Tmp files
##-----------
path_tmp = path_root+'PAHPedia/tmp/'
if not os.path.exists(path_tmp):
    os.makedirs(path_tmp)
path_conv = path_tmp+'conv/' # idlib/convolve_image.pro
if not os.path.exists(path_conv):
    os.makedirs(path_conv)
csv_ker = path_tmp+'kernelist' # idlib/conv_prog.pro

## Outputs
##---------
path_out = path_root+'PAHPedia/'+src+'/' # idlib/conv_prog.pro

## Calibrations
phot = 'IRAC4' # photometry filter
# phot = 'MIPS1' # photometry filter
path_cal = path_out+'calib/'
if not os.path.exists(path_cal):
    os.makedirs(path_cal)

## Tests
##-------
path_tests = path_root+'PAHPedia/tests/'
