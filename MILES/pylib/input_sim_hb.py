#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""

This is the interface of MILES generating following input files:

  /out/set_input.h5       -- input path
  /out/observation_MIR.h5 -- data + initial param
  /out/input_master.h5    -- general param
  /out/input_model.h5     -- model param
  /out/input_extra.h5     -- extra param

"""

import os
import math
import numpy as np

## laputan
from laputan.arrays import closest
from laputan.inout import read_fits, write_hdf5, read_hdf5

## local
from auxil import (croot, mroot,
                       res, TABLine, TABand, partuning)

## Path
##------
h5_path = mroot+'out/set_input'
dirin = croot+'../out/' # define input dir
if not os.path.exists(dirin):
    os.makedirs(dirin)
write_hdf5(h5_path, 'input dir', [dirin], verbose=True)

h5_sim = dirin+'simulation_MIR' # simulated MIR spectra
h5_obs = dirin+'observation_MIR'
h5_master = dirin+'input_master'
h5_model = dirin+'input_model'
h5_extra = dirin+'input_extra'
h5_chi2 = dirin+'fit_chi2'

program = 'fit_sim_hb'
noisy = False # verbose/debug for this routine
chi2init = True # True if using chi2 results as BB init param

##-----------------------------
##
## Write observation_MIR.h5
##
##-----------------------------
z = 0. # rest frame
wvl_inf = 5. # min wvl
wvl_sup = 20.5 # max wvl
x_inf = None # all
x_sup = None
y_inf = None
y_sup = None
# spec_unit = 'MKS' # W.m-2.Hz-1.sr-1
spec_unit = 'MJyovsr' # MJy.sr-1

## Read simulated data
wave = read_hdf5(h5_sim, 'wavelength (microns)')
data = read_hdf5(h5_sim, 'FnuOBS ('+spec_unit+')')
unc = read_hdf5(h5_sim, 'dFnuOBS ('+spec_unit+')')

## Truncate wavelength range
ind_inf = closest(wave, wvl_inf)
ind_sup = closest(wave, wvl_sup)
wave = wave[ind_inf:ind_sup]
wave = wave / (1+z)

data = data[ind_inf:ind_sup,y_inf:y_sup,x_inf:x_sup]
unc = unc[ind_inf:ind_sup,y_inf:y_sup,x_inf:x_sup]
## convert MJy/sr to W/m2/Hz/sr (MKS)
if spec_unit=='MKS':
    data = data * 1.e-20
    unc = unc * 1.e-20

## Mask NaNs
mask = np.isnan(data) * 1

## Spectroscopic modules for calibration errors
calibmod = [
            # 'IRC_NG',
            'IRS_SL2',
            'IRS_SL1',
            'IRS_LL2',
            # 'IRS_LL1',
]

## Write HDF5
##------------
write_hdf5(h5_obs, 'spectral unit', [spec_unit], verbose=True)
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
dirout = dirin
if not os.path.exists(dirout):
    os.makedirs(dirout)
verbose = 'T'
Nmcmc = 10000
NiniMC = 0 # no need for BB
calib = 'T'
robust_RMS = 'F'
robust_cal = 'F'
skew_RMS = 'F'
newseed = 'F'
dostop = 'F'
resume = 'F'
indresume = -1 # set a negative value if auto-resume
newinit = 'F'
nohi = 'F'

## Chi2 results are used as BB init param
if (chi2init):
    newinit = 'T'

## Write HDF5
##------------
write_hdf5(h5_master, 'program', [program], verbose=True)
write_hdf5(h5_master, 'output dir', [dirout], verbose=noisy)
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
write_hdf5(h5_master, 'resume', [resume], append=True, verbose=noisy)
write_hdf5(h5_master, 'indresume', [indresume], append=True, verbose=noisy)
write_hdf5(h5_master, 'nohi', [nohi], append=True, verbose=noisy)

##--------------------------------
##
## Write input_model.h5
##
##--------------------------------
labQ = ['BE_Z96               ', # 10
        'BE_Z96               ', # 10
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
        'Main 6.2 (1) ', # 7*
        'Main 6.2 (2) ', # 8
        'Plateau 7.7  ', # 12
        'Main 7.7 (1) ', # 13*
        'Main 7.7 (2) ', # 14
        'Main 8.6     ', # 16*
        'Small 11.0   ', # 19
        'Main 11.2    ', # 20*
        'Plateau 11.3 ', # 21
        'Main 12.7 (1)', # 24
        'Main 12.7 (2)', # 25
        'Plateau 17.0 '] # 30
labE = ['D03']

refB = ['Main 11.2    ']
refw = 15.0 # ISO/LW9 (15 µm) pure continuum emission

ALline = False
ALband = True

if (ALline):
    labL = []
    for tabL in TABLine:
        if (tabL['wave']>wave[0] and tabL['wave']<wave[-1]):
            labL.append(tabL['label'])
if (ALband):
    labB = []
    for tabB in TABand:
        if (tabB['wave']>wave[0] and tabB['wave']<wave[-1]):
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
            ## lnFcont, lnRband(ref), lnFstar,
            ##----------------------------------------------
            dict([ ('namall','lnFcont'),('fixed','F'),('hyper','T'),]),
            
            dict([ ('namall','lnRline'),('fixed','F'),('hyper','T'),]),
            
            dict([ ('namall','lnRband'),('fixed','F'),('hyper','T'),]),
            
            dict([ ('namall','lnFstar'),('fixed','F'),('hyper','T'),]),
            
            ## Intensive param:
            ## lnT, lnRline, lnRband, lnAv,
            ## (fixed) Cline, Cband, Wline, WSband, WLband,
            ##----------------------------------------------
            dict([ ('namall','lnT'),
                   ('fixed','F'),('hyper','T'),
                   ('limited',('T','T')),
                   ('limits',(np.log(50.),np.log(500.))),
                   ('hyper','T'),
            ]), # LOG( (50,500) K )
            
            # dict([ ('namall','Cline'),('fixed','F'),('hyper','T'),]),
            
            ## Main 6.2 (1)
            dict([ ('name','Cband'+str(labB.index('Main 6.2 (1)')+1)),
                   ('fixed','F'),('hyper','T') ]),
            dict([ ('name','WSband'+str(labB.index('Main 6.2 (1)')+1)),
                   ('fixed','F'),('hyper','T'), ]),
            dict([ ('name','WLband'+str(labB.index('Main 6.2 (1)')+1)),
                   ('fixed','F'),('hyper','T') ]),
            ## Main 7.7 (1)
            dict([ ('name','Cband'+str(labB.index('Main 7.7 (1)')+1)),
                   ('fixed','F'),('hyper','T') ]),
            dict([ ('name','WSband'+str(labB.index('Main 7.7 (1)')+1)),
                   ('fixed','F'),('hyper','T') ]),
            dict([ ('name','WLband'+str(labB.index('Main 7.7 (1)')+1)),
                   ('fixed','F'),('hyper','T') ]),
            ## Main 8.6
            dict([ ('name','Cband'+str(labB.index('Main 8.6')+1)),
                   ('fixed','F'),('hyper','T') ]),
            dict([ ('name','WSband'+str(labB.index('Main 8.6')+1)),
                   ('fixed','F'),('hyper','T') ]),
            dict([ ('name','WLband'+str(labB.index('Main 8.6')+1)),
                   ('fixed','F'),('hyper','T') ]),
            ## Main 11.2
            dict([ ('name','Cband'+str(labB.index('Main 11.2')+1)),
                   ('fixed','F'),('hyper','T') ]),
            dict([ ('name','WSband'+str(labB.index('Main 11.2')+1)),
                   ('fixed','F'),('hyper','T') ]),
            dict([ ('name','WLband'+str(labB.index('Main 11.2')+1)),
                   ('fixed','F'),('hyper','T') ]),
            
            dict([ ('namall','lnAv'),('fixed','F'),('hyper','T'), ]),
            # dict([ ('namall','lnAv'),('fixed','T'),('value','0.5'), ]), # LOG( exp(.5) mag )
            
            dict() ]

## parinfo
##---------
Ncont = len(labQ)
Nline = len(labL)
Nband = len(labB)
Nextc = len(labE)
Nstar = 1
Npar = 2*Ncont + 3*Nline + 4*Nband + Nextc + Nstar

## Set indL (index for labL)
indL = []
for lab in labL:
    for tabL in TABLine:
        if tabL['label']==lab.rstrip():
            indL.append(TABLine.index(tabL))
## Set indB (index for labB)
indB = []
for lab in labB:
    for tabB in TABand:
        if tabB['label']==lab.rstrip():
            indB.append(TABand.index(tabB))

## Default setting
name = np.empty((Npar,), dtype=('<U30')) # str length <30
namall = np.empty((Npar,), dtype=('<U30')) # specific for partuning
comp = np.empty((Npar,), dtype=('<U30'))
fixed = np.array(['T' for i in range(Npar)]) # fixed by default
limited = np.array([('F','F') for i in range(Npar)])
limits = np.array([(0.,0.) for i in range(Npar)])
model = np.array(['T' for i in range(Npar)])
hyper = np.array(['F' for i in range(Npar)])
tied = np.empty((Npar,), dtype=('<U30'))
value = np.array([0. for i in range(Npar)])

## Param assignment
i0 = 0
for i in range(Ncont):
    name[i0+2*i] = 'lnFcont'+str(i+1)
    namall[i0+2*i] = 'lnFcont'
    value[i0+2*i] = 0. # 1 [W/m2/sr]
    name[i0+2*i+1] = 'lnT'+str(i+1)
    namall[i0+2*i+1] = 'lnT'
    value[i0+2*i+1] = 4. # 54.60 [K]
for i in range(2*Ncont):
    comp[i0+i] = 'CONT'
i0 += 2*Ncont
for i in range(Nline):
    name[i0+3*i] = 'lnRline'+str(i+1)
    namall[i0+3*i] = 'lnRline'
    value[i0+3*i] = 0.
    name[i0+3*i+1] = 'Cline'+str(i+1)
    namall[i0+3*i+1] = 'Cline'
    value[i0+3*i+1] = TABLine[indL[i]]['wave']
    name[i0+3*i+2] = 'Wline'+str(i+1)
    namall[i0+3*i+2] = 'Wline'
    for r in res:
        if value[i0+3*i+1]<5.:
            if r['name']=='AKARI_Ns':
                value[i0+3*i+2] = r['dwovw'] * value[i0+3*i+1]
        elif value[i0+3*i+1]<12.:
            if r['name']=='SL':
                value[i0+3*i+2] = r['dwovw'] * value[i0+3*i+1]
        else:
            if r['name']=='LL':
                value[i0+3*i+2] = r['dwovw'] * value[i0+3*i+1]
for i in range(3*Nline):
    comp[i0+i] = 'LINE'
i0 += 3*Nline
for i in range(Nband):
    name[i0+4*i] = 'lnRband'+str(i+1)
    namall[i0+4*i] = 'lnRband'
    value[i0+4*i] = 0.
    name[i0+4*i+1] = 'Cband'+str(i+1)
    namall[i0+4*i+1] = 'Cband'
    value[i0+4*i+1] = TABand[indB[i]]['wave']
    name[i0+4*i+2] = 'WSband'+str(i+1)
    namall[i0+4*i+2] = 'WSband'
    value[i0+4*i+2] = TABand[indB[i]]['sigmaS']
    name[i0+4*i+3] = 'WLband'+str(i+1)
    namall[i0+4*i+3] = 'WLband'
    value[i0+4*i+3] = TABand[indB[i]]['sigmaL']
for i in range(4*Nband):
    comp[i0+i] = 'BAND'
i0 += 4*Nband
for i in range(Nextc):
    name[i0+i] = 'lnAv'+str(i+1)
    namall[i0+i] = 'lnAv'
    value[i0+i] = 0. # 1 [mag]
    comp[i0+i] = 'EXTC'
i0 += Nextc
for i in range(Nstar):
    name[i0+i] = 'lnFstar'+str(i+1)
    namall[i0+i] = 'lnFstar'
    value[i0+i] = 0. # 1 [W/m2]
    comp[i0+i] = 'STAR'
    
## Param tuning
partuning(dictune, Ncont, Nline, Nband, Nextc,
          name, fixed, limited, limits, model, hyper, tied, value)
    
## Write HDF5
##------------
write_hdf5(h5_model, 'program', [program], verbose=True)
write_hdf5(h5_model, 'label cont', labQ, append=True, verbose=noisy)
write_hdf5(h5_model, 'label band', labB, append=True, verbose=noisy)
write_hdf5(h5_model, 'label line', labL, append=True, verbose=noisy)
write_hdf5(h5_model, 'label extc', labE, append=True, verbose=noisy)
write_hdf5(h5_model, 'ref band', refB, append=True, verbose=noisy)
write_hdf5(h5_model, 'ref wavelength', refw, append=True, verbose=noisy)
write_hdf5(h5_model, 'parinfo name', name, append=True, verbose=noisy)
write_hdf5(h5_model, 'parinfo comp', comp, append=True, verbose=noisy)
write_hdf5(h5_model, 'parinfo fixed', fixed, append=True, verbose=noisy)
write_hdf5(h5_model, 'parinfo limited', limited, append=True, verbose=noisy)
write_hdf5(h5_model, 'parinfo limits', limits, append=True, verbose=noisy)
write_hdf5(h5_model, 'parinfo model', model, append=True, verbose=noisy)
write_hdf5(h5_model, 'parinfo hyper', hyper, append=True, verbose=noisy)
write_hdf5(h5_model, 'parinfo tied', tied, append=True, verbose=noisy)
write_hdf5(h5_model, 'parinfo value', value, append=True, verbose=noisy)

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
write_hdf5(h5_obs, 'Initial parameter label', name, append=True, verbose=noisy)
val2 = np.repeat(value[:,np.newaxis], data.shape[1], axis=1) # expand Ny
val3 = np.repeat(val2[:,:,np.newaxis], data.shape[2], axis=2) # expand Nx
write_hdf5(h5_obs, 'Initial parameter value', val3, append=True, verbose=noisy)
if (chi2init):
    write_hdf5(h5_obs, 'Chi2init', ['T'], append=True, verbose=True)
    parini = read_hdf5(h5_chi2, 'Best fitted parameter value')
    write_hdf5(h5_obs, 'Chi2 fitted parameter value', parini,
               append=True, verbose=True)
else:
    write_hdf5(h5_obs, 'Chi2init', ['F'], append=True, verbose=True)
