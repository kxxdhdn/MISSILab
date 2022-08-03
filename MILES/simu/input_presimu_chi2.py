#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""

This is the interface of MILES generating following files:

  /out/set_input.h5       -- input path
  /out/observation_MIR.h5 -- data + initial param
  /out/input_master.h5    -- general param
  /out/input_model.h5     -- model param
  /out/input_extra.h5     -- extra param

"""

import os
import math
import numpy as np

## rapyuta
from rapyuta.arrays import closest, pix2sup
from rapyuta.inout import read_fits, write_hdf5, read_hdf5

## local
from auxil import (croot, mroot, TABLine, TABand, partuning)
         
## Path
##------
h5_path = mroot+'out/set_input'
dirin = croot+'../out/' # define input dir
os.makedirs(dirin, exist_ok=True)
write_hdf5(h5_path, 'input dir', [dirin], verbose=True)

h5_obs = dirin+'observation_MIR'
h5_master = dirin+'input_master'
h5_model = dirin+'input_model'
h5_extra = dirin+'input_extra'

program = 'presimu_chi2'
noisy = False # verbose/debug for this routine

##-----------------------------
##
## Write observation_MIR.h5
##
##-----------------------------

## Read FITS
##-----------
spec_name = [['A1']]
xscale, yscale = 3, 6
fits_obs = mroot+'lib/M82_A1_7'
fits_unc = fits_obs+'_unc' # unc
ds = read_fits(fits_obs, fits_unc)
data0 = ds.data[:,0,0][:,np.newaxis,np.newaxis]
unc0 = ds.unc[:,0,0][:,np.newaxis,np.newaxis]
wave0 = ds.wave

## Data
##------
z = 0.00068
winf = 2.50 # min wvl
wsup = 20.00 # max wvl
# spec_unit = 'MKS' # W.m-2.Hz-1.sr-1
spec_unit = 'MJyovsr' # MJy.sr-1

## Truncate wavelength range
ind_inf = closest(wave0, winf)
ind_sup = closest(wave0, wsup)
data = data0[ind_inf:ind_sup]
unc = unc0[ind_inf:ind_sup]
wave = wave0[ind_inf:ind_sup]
wave = wave / (1+z)

## convert MJy/sr to W/m2/Hz/sr (MKS)
if spec_unit=='MKS':
    data = data * 1.e-20
    unc = unc * 1.e-20

## Mask NaNs
mask = np.isnan(data) * 1
## Mask edge
i1 = closest(wave0, 5, 'left')
i2 = closest(wave0, 5.35, 'right')
mask[i1:i2] = 1

## Spectroscopic modules for calibration errors
calibmod = [
            'IRS_SL2', # ref
            'IRC_NG',
            'IRS_SL1',
            'IRS_LL2',
            # 'IRS_LL1',
]

## Write HDF5
##------------
write_hdf5(h5_obs, 'spectral unit', [spec_unit], verbose=True)
write_hdf5(h5_obs, 'spectrum labels', spec_name, append=True, verbose=noisy)
write_hdf5(h5_obs, 'wavelength (microns)', wave, append=True, verbose=noisy)
write_hdf5(h5_obs, 'FnuOBS ('+spec_unit+')', data, append=True, verbose=noisy)
write_hdf5(h5_obs, 'dFnuOBS ('+spec_unit+')', unc, append=True, verbose=noisy)
write_hdf5(h5_obs, 'NaN mask', mask, append=True, verbose=noisy)
write_hdf5(h5_obs, 'spectroscopic module labels', calibmod, append=True, verbose=noisy)

##--------------------------------
##
## Write input_master.h5
##
##--------------------------------
dirout = dirin # define output dir
os.makedirs(dirout, exist_ok=True)
verbose = 'T'
Nmcmc = 10000 # no need for chi2
NiniMC = 10
calib = 'F'
robust_RMS = 'F'
robust_cal = 'F'
skew_RMS = 'F'
newseed = 'F'
dostop = 'F'
newinit = 'F'

## Write HDF5
##------------
write_hdf5(h5_master, 'program', [program], verbose=True)
write_hdf5(h5_master, 'output dir', [dirout], append=True, verbose=noisy)
write_hdf5(h5_master, 'spectral unit', [spec_unit], append=True, verbose=noisy)
write_hdf5(h5_master, 'verbose', [verbose], append=True, verbose=noisy)
write_hdf5(h5_master, 'Nmcmc', [Nmcmc], append=True, verbose=noisy)
write_hdf5(h5_master, 'NiniMC', [NiniMC], append=True, verbose=noisy)
write_hdf5(h5_master, 'calib', [calib], append=True, verbose=noisy)
write_hdf5(h5_master, 'robust_RMS', [robust_RMS], append=True, verbose=noisy)
write_hdf5(h5_master, 'robust_cal', [robust_cal], append=True, verbose=noisy)
write_hdf5(h5_master, 'skew_RMS', [skew_RMS], append=True, verbose=noisy)
write_hdf5(h5_master, 'newseed', [newseed], append=True, verbose=noisy)
write_hdf5(h5_master, 'newinit', [newinit], append=True, verbose=noisy)
write_hdf5(h5_master, 'dostop', [dostop], append=True, verbose=noisy)

##--------------------------------
##
## Write input_model.h5
##
##--------------------------------
labQ = [
        'BE_Z96               ', # 10
        # 'BE_Z96               ', # 10
        # 'BE_Z96               ', # 10
        'Gra_D03              ', # 12
        # 'Gra_D03              ', # 12
        'Sil_D03              ', # 23
        ]
labL = ['Bra   ', # 1
        # 'H2S7  ', # 3
        # 'Huc   ', # 5
        # 'H2S5  ', # 9
        'ArII  ', # 11
        'ArIII1', # 20
        'H2S3  ', # 23
        # 'HeII2 ', # 24
        'SIV   ', # 26
        'H2S2  ', # 27
        'NeII  ', # 29
        # 'NeV1  ', # 32
        'NeIII1', # 33
        'H2S1  ', # 34
        'SIII1 '] # 36
        # 'ArIII2', # 38
        # 'FeII1 ', # 41
        # 'H2S0. ', # 42
        # 'SIII2 ', # 43
        # 'NeIII2'] # 46
labB = ['Main 3.3     ', # 1
        'Main 3.4     ', # 2
        'Small 3.5    ', # 3
        'Main 6.2 (1) ', # 8*
        'Main 6.2 (2) ', # 9
        'Plateau 7.7  ', # 13
        'Main 7.7 (1) ', # 14*
        'Main 7.7 (2) ', # 15
        'Main 8.6     ', # 17*
        'Small 11.0   ', # 20
        'Main 11.2    ', # 21*
        'Plateau 11.3 ', # 22
        'Main 12.7 (1)', # 25
        'Main 12.7 (2)', # 26
        'Plateau 17.0 '] # 31
labE = ['D03']
labS = ['BB1']

refB = ['Main 11.2    ']
refw = 15.0 # ISO/LW9 (15 µm) pure continuum emission

ALline = False
ALband = False

if (ALline):
    labL = []
    for tabL in TABLine:
        if (tabL['wave']>wave[0] and tabL['wave']<wave[-1]):
            labL.append(tabL['label'])
if (ALband):
    labB = []
    for tabB in TABand:
        if (tabB['wave']>wave[0] and tabB['wave']<wave[-1]
            and tabB['label']!='Small 5.2'): # exclude 5.2 um band
            labB.append(tabB['label'])

## dictune, as input of the partuning function,
## is a list of dict containing param tuning info
dictune = [ dict([ ('name','default'),
                   ('namall','default'),
                   ('fixed','T'),
                   ('limits',(0.,0.)),
                   ('limited',('F','F')),
                   ('model','T'),
                   ('hyper','F'),
                   ('tied',''),
                   ('value',0.), ]),
            
            ## Add your config below
            ##=======================

            ## Extensive param:
            ## lnFcont, lnRband(ref), lnLstar,
            ##----------------------------------------------
            dict([ ('namall','lnFcont'),('fixed','F'),]),
            
            dict([ ('namall','lnRline'),('fixed','F'),]),
            
            dict([ ('namall','lnRband'),('fixed','F'),]),
            
            dict([ ('namall','lnLstar'),('fixed','F'),]),
            
            ## Intensive param:
            ## lnT, lnRline, lnRband, lnAv,
            ## (fixed) Cline, Cband, Wline, WSband, WLband,
            ##----------------------------------------------
            dict([ ('namall','lnT'),
                   ('fixed','F'),
                   ('limited',('T','T')),
                   ('limits',(np.log(50.),np.log(500.))),
            ]), # LOG( (50,500) K )
            
            # dict([ ('namall','Cline'),('fixed','F'),]),

            # ## Main 3.3
            # dict([ ('name','Cband'+str(labB.index('Main 3.3     ')+1)),('fixed','F') ]),
            # dict([ ('name','WSband'+str(labB.index('Main 3.3     ')+1)),('fixed','F') ]),
            # dict([ ('name','WLband'+str(labB.index('Main 3.3     ')+1)),('fixed','F') ]),
            # ## Main 3.4
            # dict([ ('name','Cband'+str(labB.index('Main 3.4     ')+1)),('fixed','F') ]),
            # dict([ ('name','WSband'+str(labB.index('Main 3.4     ')+1)),('fixed','F') ]),
            # dict([ ('name','WLband'+str(labB.index('Main 3.4     ')+1)),('fixed','F') ]),
            # ## Small 3.5
            # dict([ ('name','Cband'+str(labB.index('Small 3.5    ')+1)),('fixed','F') ]),
            # dict([ ('name','WSband'+str(labB.index('Small 3.5    ')+1)),('fixed','F') ]),
            # dict([ ('name','WLband'+str(labB.index('Small 3.5    ')+1)),('fixed','F') ]),
            # ## Main 6.2 (1)
            # dict([ ('name','Cband'+str(labB.index('Main 6.2 (1) ')+1)),('fixed','F') ]),
            # ## Main 6.2 (2)
            # dict([ ('name','Cband'+str(labB.index('Main 6.2 (2) ')+1)),('fixed','F') ]),
            # ## Plateau 7.7
            # dict([ ('name','Cband'+str(labB.index('Plateau 7.7  ')+1)),('fixed','F') ]),
            # ## Main 7.7 (1)
            # dict([ ('name','Cband'+str(labB.index('Main 7.7 (1) ')+1)),('fixed','F') ]),
            # ## Main 7.7 (2)
            # dict([ ('name','Cband'+str(labB.index('Main 7.7 (2) ')+1)),('fixed','F') ]),
            # ## Main 8.6
            # dict([ ('name','Cband'+str(labB.index('Main 8.6     ')+1)),('fixed','F') ]),
            # ## Small 11.0
            # dict([ ('name','Cband'+str(labB.index('Small 11.0   ')+1)),('fixed','F') ]),
            # ## Main 11.2
            # dict([ ('name','Cband'+str(labB.index('Main 11.2    ')+1)),('fixed','F') ]),
            # ## Plateau 11.3
            # dict([ ('name','Cband'+str(labB.index('Plateau 11.3 ')+1)),('fixed','F') ]),
            # ## Main 12.7 (1)
            # dict([ ('name','Cband'+str(labB.index('Main 12.7 (1)')+1)),('fixed','F') ]),
            # ## Main 12.7 (2)
            # dict([ ('name','Cband'+str(labB.index('Main 12.7 (2)')+1)),('fixed','F') ]),
            # ## Plateau 17.0
            # dict([ ('name','Cband'+str(labB.index('Plateau 17.0 ')+1)),('fixed','F') ]),
            
            dict([ ('namall','lnAv'),('fixed','F'),]),
            
            dict() ]

## Parameter tuning
ds = partuning(dictune, labQ, labL, labB, labE, labS)

## Write HDF5
##------------
write_hdf5(h5_model, 'program', [program], verbose=True)
write_hdf5(h5_model, 'label cont', labQ, append=True, verbose=noisy)
write_hdf5(h5_model, 'label band', labB, append=True, verbose=noisy)
write_hdf5(h5_model, 'label line', labL, append=True, verbose=noisy)
write_hdf5(h5_model, 'label extc', labE, append=True, verbose=noisy)
write_hdf5(h5_model, 'label star', labS, append=True, verbose=noisy)
write_hdf5(h5_model, 'ref band', refB, append=True, verbose=noisy)
write_hdf5(h5_model, 'ref wavelength', refw, append=True, verbose=noisy)
write_hdf5(h5_model, 'parinfo name', ds.name, append=True, verbose=noisy)
write_hdf5(h5_model, 'parinfo comp', ds.comp, append=True, verbose=noisy)
write_hdf5(h5_model, 'parinfo fixed', ds.fixed, append=True, verbose=noisy)
write_hdf5(h5_model, 'parinfo limited', ds.limited, append=True, verbose=noisy)
write_hdf5(h5_model, 'parinfo limits', ds.limits, append=True, verbose=noisy)
write_hdf5(h5_model, 'parinfo model', ds.model, append=True, verbose=noisy)
write_hdf5(h5_model, 'parinfo hyper', ds.hyper, append=True, verbose=noisy)
write_hdf5(h5_model, 'parinfo tied', ds.tied, append=True, verbose=noisy)
write_hdf5(h5_model, 'parinfo value', ds.value, append=True, verbose=noisy)

##--------------------------------
##
## Write input_extra.h5
##
##--------------------------------
Nextra = 0

## parextinfo
##------------

## Write HDF5
##------------
write_hdf5(h5_extra, 'program', [program], verbose=True)
write_hdf5(h5_extra, 'Nextra', [Nextra], append=True, verbose=noisy)


##-----------------------------
##
## Append observation_MIR.h5
##
##-----------------------------

## These init param are supposed to be the default param in the fitting model,
## with the possibility of reasonable modifications by this script.
write_hdf5(h5_obs, 'Initial parameter label', ds.name, append=True, verbose=True)
val2 = np.repeat(ds.value[:,np.newaxis], data.shape[1], axis=1) # expand Ny
val3 = np.repeat(val2[:,:,np.newaxis], data.shape[2], axis=2) # expand Nx
write_hdf5(h5_obs, 'Initial parameter value', val3, append=True, verbose=True)
write_hdf5(h5_obs, 'Chi2init', ['F'], append=True, verbose=True)
