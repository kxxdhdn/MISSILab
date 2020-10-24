#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""

This is the interface of the MILES input files as follows:

  observations_fitMIR.h5   -- data
  input_fitMIR_master.h5   -- general param
  input_fitMIR_model.h5    -- model param
  input_fitMIR_extra.h5    -- extra param
  input_fitMIR_analysis.h5 -- param for analysis

"""

import os
import math
import numpy as np

## astylo
from astylo.arrlib import closest
from astylo.iolib import read_fits, write_hdf5, read_hdf5

## local
from utilities import TABLine, TABand, partuning
         
## Path
##------
mroot = os.path.dirname(os.path.abspath(__file__))+'/..'
path_in = mroot+'/tests'
path_data = path_in+'/dat'
h5_obs = path_data+'/observations_fitMIR'
h5_master = path_in+'/input_fitMIR_master'
h5_model = path_in+'/input_fitMIR_model'
h5_extra = path_in+'/input_fitMIR_extra'
h5_analysis = path_in+'/input_fitMIR_analysis'

##-----------------------------
## Set parameters
##-----------------------------
noisy = False # verbose/debug for this routine

##-----------------------------
## Write observations_fitMIR
##-----------------------------
fits_obs = path_data+'/M83' # obs
fits_unc = fits_obs+'_unc' # unc
wvl_inf = 1. # min wvl
wvl_sup = 21. # max wvl

## Read FITS
dset = read_fits(fits_obs, fits_unc)

## Truncate wavelength range
ind_inf = closest(dset.wave, wvl_inf)
ind_sup = closest(dset.wave, wvl_sup)
wave = dset.wave[ind_inf:ind_sup]
data = dset.data[ind_inf:ind_sup,:,:]
unc = dset.unc[ind_inf:ind_sup,:,:]
## Mask NaNs
mask = np.isnan(dset.data) * 1

## Write HDF5
##------------
write_hdf5(h5_obs, 'Wavelength (microns)', wave, verbose=noisy)
write_hdf5(h5_obs, 'FnuOBS (MJyovsr)', data, append=True, verbose=noisy)
write_hdf5(h5_obs, 'dFnuOBS (MJyovsr)', unc, append=True, verbose=noisy)
write_hdf5(h5_obs, 'NaN mask', mask, append=True, verbose=noisy)

##--------------------------------
## Write input_fitMIR_master.h5
##--------------------------------
verbose = 'T'
Nmcmc = 200
NiniMC = 100
calib = 'F'
robust_RMS = 'F'
robust_cal = 'F'
skew_RMS = 'F'
newseed = 'F'
newinit = 'F'
dostop = 'F'

## Write HDF5
##------------
write_hdf5(h5_master, 'verbose', [verbose], verbose=noisy)
write_hdf5(h5_master, 'Nmcmc', [Nmcmc], append=True, verbose=noisy)
write_hdf5(h5_master, 'NiniMC', [NiniMC], append=True, verbose=noisy)
write_hdf5(h5_master, 'calib', [calib], append=True, verbose=noisy)
write_hdf5(h5_master, 'robust_RMS', [robust_RMS], append=True, verbose=noisy)
write_hdf5(h5_master, 'robust_cal', [robust_cal], append=True, verbose=noisy)
write_hdf5(h5_master, 'skew_RMS', [skew_RMS], append=True, verbose=noisy)
write_hdf5(h5_master, 'newseed', [newseed], append=True, verbose=noisy)
write_hdf5(h5_master, 'newinit', [newinit], append=True, verbose=noisy)
write_hdf5(h5_master, 'dostop', [dostop], append=True, verbose=noisy)
# write_hdf5(h5_master, '', , append=True, verbose=noisy)

##--------------------------------
## Write input_fitMIR_model.h5
##--------------------------------
spec_unit = 'MJy.sr-1'

labQ = ['ACH2_Z96             ', # 9
        'BE_Z96               ', # 10
        'Sil_D03              '] # 23
labL = ['H2S7  ', # 3
        'H2S5  ', # 9
        'ArII  ', # 11
        'ArIII1', # 20
        'H2S3  ', # 23
        'HeII2 ', # 24
        'SIV   ', # 26
        'H2S2  ', # 27
        'NeII  ', # 29
        'NeIII1', # 33
        'H2S1  ', # 34
        'SIII1 '] # 36
labB = ['Main 3.3     ', # 1
        'Main 3.4     ', # 2
        'Main 6.2 (1) ', # 7
        'Main 6.2 (2) ', # 8
        'Plateau 7.7  ', # 12
        'Main 7.7 (1) ', # 13
        'Main 7.7 (2) ', # 14
        'Main 8.6     ', # 16
        'Small 11.0   ', # 19
        'Main 11.2    ', # 20
        'Plateau 11.3 ', # 21
        'Main 12.7 (1)', # 24
        'Main 12.7 (2)', # 25
        'Plateau 17.0 '] # 30

ALline = False
ALband = False

## dictune, as input of the partuning function,
## is a list of dict containing param tuning info
dictune = [ dict([ ('name','default'),
                   ('namall','default'),
                   ('fixed','T'),
                   ('limited',('F','F')),
                   ('limits',(0.,0.)),
                   ('hyper','F'),
                   ('tied',''),
                   ('value',0.) ]),
            ## Add your config below
            dict([ ('namall','lnMovd2'),
                   ('fixed','F'),
                   ('limited',('T','F')) ]),
            dict([ ('namall','lnT'),
                   ('fixed','F'),
                   ('limited',('T','T')),
                   ('limits',(np.log(50.),np.log(500.))) ]), # LOG( (50,500) K ) 
            dict([ ('namall','lnIline'),
                   ('fixed','F'),
                   ('limited',('T','T')) ]),
            dict([ ('namall','Cline'),
                   ('fixed','F'),
                   ('limited',('T','T')) ]),
            dict([ ('namall','lnIband'),
                   ('fixed','F'),
                   ('limited',('T','T')) ]),
            dict([ ('namall','Cband'),
                   ('fixed','F'),
                   ('limited',('T','F')) ]),
            dict([ ('name','lnFstar'),
                   ('fixed','F'),
                   ('limited',('T','F')) ]),
            dict() ]

## parinfo
##---------
Ncont = len(labQ)
if (ALline):
    labL = []
    for tabL in TABLine:
        labL.append(tabL['label'])
Nline = len(labL)
if (ALband):
    labB = []
    for tabB in TABand:
        labB.append(tabB['label'])
Nband = len(labB)
Npar = 2*Ncont + 3*Nline + 4*Nband + 2

## Default setting
name = np.empty((Npar,), dtype=('<U30')) # str length <30
namall = np.empty((Npar,), dtype=('<U30')) # specific for partuning
comp = np.empty((Npar,), dtype=('<U30'))
fixed = np.array(['T' for i in range(Npar)])
limited = np.array([('F','F') for i in range(Npar)])
limits = np.array([(0.,0.) for i in range(Npar)])
hyper = np.array(['F' for i in range(Npar)])
tied = np.empty((Npar,), dtype=('<U30'))
value = np.array([0. for i in range(Npar)])

## Param assignment
i0 = 0
for i in range(Ncont):
    name[i0+2*i] = 'lnMovd2'+str(i+1)
    namall[i0+2*i] = 'lnMovd2'
    name[i0+2*i+1] = 'lnT'+str(i+1)
    namall[i0+2*i+1] = 'lnT'
for i in range(2*Ncont):
    comp[i0+i] = 'CONT'
i0 += 2*Ncont
for i in range(Nline):
    name[i0+3*i] = 'lnIline'+str(i+1)
    namall[i0+3*i] = 'lnIline'
    name[i0+3*i+1] = 'Cline'+str(i+1)
    namall[i0+3*i+1] = 'Cline'
    name[i0+3*i+2] = 'Wline'+str(i+1)
    namall[i0+3*i+2] = 'Wline'
for i in range(3*Nline):
    comp[i0+i] = 'LINE'
i0 += 3*Nline
for i in range(Nband):
    name[i0+4*i] = 'lnIband'+str(i+1)
    namall[i0+4*i] = 'lnIband'
    name[i0+4*i+1] = 'Cband'+str(i+1)
    namall[i0+4*i+1] = 'Cband'
    name[i0+4*i+2] = 'WSband'+str(i+1)
    namall[i0+4*i+2] = 'WSband'
    name[i0+4*i+3] = 'WLband'+str(i+1)
    namall[i0+4*i+3] = 'WLband'
for i in range(4*Nband):
    comp[i0+i] = 'BAND'
i0 += 4*Nband
name[i0] = 'lnAv'
namall[i0] = 'lnAv'
comp[i0] = 'PABS'
i0 += 1
name[i0] = 'lnFstar'
namall[i0] = 'lnFstar'
comp[i0] = 'STAR'

## Param tuning
partuning(dictune, Ncont, Nline, Nband,
          name, fixed, limited, limits, hyper, tied, value)

## Write HDF5
##------------
write_hdf5(h5_model, 'spec_unit', [spec_unit], verbose=noisy)
write_hdf5(h5_model, 'labQ', labQ, append=True, verbose=noisy)
write_hdf5(h5_model, 'labB', labB, append=True, verbose=noisy)
write_hdf5(h5_model, 'labL', labL, append=True, verbose=noisy)
write_hdf5(h5_model, 'parinfo name', name, append=True, verbose=noisy)
write_hdf5(h5_model, 'parinfo comp', comp, append=True, verbose=noisy)
write_hdf5(h5_model, 'parinfo fixed', fixed, append=True, verbose=noisy)
write_hdf5(h5_model, 'parinfo limited', limited, append=True, verbose=noisy)
write_hdf5(h5_model, 'parinfo limits', limits, append=True, verbose=noisy)
write_hdf5(h5_model, 'parinfo hyper', hyper, append=True, verbose=noisy)
write_hdf5(h5_model, 'parinfo tied', tied, append=True, verbose=noisy)
write_hdf5(h5_model, 'parinfo value', value, append=True, verbose=noisy)

##--------------------------------
## Write input_fitMIR_extra.h5
##--------------------------------
Nextra = 0

## parextinfo
##------------

## Write HDF5
##------------
write_hdf5(h5_extra, 'Nextra', [Nextra], verbose=noisy)

##--------------------------------
## Write input_fitMIR_analysis.h5
##--------------------------------

## Write HDF5
##------------

