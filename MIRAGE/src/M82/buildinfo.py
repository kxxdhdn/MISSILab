#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""

SOURCE: M82 (NGC 3034)

"""

import os


##----------------------------------------------------------

##                    Inputs & Outputs

##----------------------------------------------------------
src = 'M82'
Nmc = 100
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
parobs = [
    ## [Name,obsid,,,,,,spec,,imref,,,,,,,,,,,slit,Ny,Nsub,Nx,colors,markers]
    ['A1_7','3390001.2','NG','F011100338_N002','Nh',28,7,3,'r','*'], # A1-7
    ['B1_7','3390003.1','NG','F011100379_N002','Nh',28,7,3,'orange','*'], # B1-7
    ['C1_7','3390002.1','NG','F011100795_N002','Nh',28,7,3,'y','*'], # C1-7
    ['D1_4','5125405.1','NG','F010116950_N002','Nh',24,4,3,'g','d'], # D1-4
    ['E1_4','5125401.1','NG','F010117172_N002','Ns',24,4,5,'lime','d'], # E1-4
    ['F1_2','5125405.1','NG','F010116950_N002','Ns',24,2,5,'pink','v'], # F1-2
    ['G1_2','3390003.1','NG','F011100379_N002','Ns',24,2,5,'m','v'], # G1-2
    ['H1_2','3390001.2','NG','F011100338_N002','Ns',24,2,5,'c','v'], # H1-2
    ['I1_2','3390002.1','NG','F011100795_N002','Ns',24,2,5,'b','v'], # I1-2
    ['J1_2','5125401.1','NG','F010117172_N002','Nh',24,2,3,'b','^'], # J1-2
    ['K1_2','5125403.1','NG','F010116924_N002','Ns',24,2,5,'m','^'], # K1-2
    ['L1_2','5125403.1','NG','F010116924_N002','Nh',24,2,3,'pink','^'], # L1-2
    ['M1_2','5124077.1','NG','F007174142_N002','Ns',24,2,5,'c','^'], # M1-2
    ['N1_2','5124077.1','NG','F007174142_N002','Nh',24,2,3,'c','s'], # N1-2
    # ['3390001.1','NG','F011100297_N002','Ns',24,2,5], #
    # ['3390001.1','NG','F011100297_N002','Nh',24,7,3], #
    # ['3390002.2','NG','F011176073_N002'], # Matching failed
    # ['5125402.1','NP','F010117172_N002'], # Matching failed
    # ['5125404.1','NP','F010117338_N002'], # Matching failed
    # ['5125406.1','NP','F010117086_N002'], # NP has 68 wvl instead of 259
]

fits_irc = []
out_irc = []
out_irs = []
colors = ['k']
markers = []
for obs in parobs:
    fits_irc.append(path_build+obs[0])
    out_irc.append(path_out+src+'_'+obs[0]+'_IRC')
    out_irs.append(path_out+src+'_'+obs[0]+'_IRS')
    colors.append(obs[8])
    markers.append(obs[9])

filog = path_out+'build_history_'


##----------------------------------------------------------

##                  IRS data (via CUBISM)

##----------------------------------------------------------
# out_irs = path_out+src+'_IRS'

# chnl = ['SH', 'LH']
chnl = ['SL2', 'SL1', 'LL2', 'LL1', 'SL3', 'LL3']

# sub_SL = ['04','06S','06N','08','08c','09N3','09N2']
# sub_LL = ['04','05','06','08','09N3','09N5','09N2']
sub_SL = ['04','06S','06N','08c','09N3']
sub_LL = ['04','05','06','08','09N3']
sub_SH = []
sub_LH = []

##----------------------------------------------------------

##                   Data homogenisation

##----------------------------------------------------------
# coadd_tool = 'swarp'
coadd_tool = 'reproject'

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
