#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""

This is the visualization of presimulation.h5 (chi2 fit)

"""

import os, pathlib
import numpy as np
import matplotlib.pyplot as plt

## rapyuta
from rapyuta.arrays import closest
from rapyuta.inout import read_hdf5
from rapyuta.plots import pplot

## local
from auxil import croot
                       
## Path
##------
path_out = croot+'/../out/'
filout = path_out+'presimulation'
path_fig = path_out+'Figures/'
if not os.path.exists(path_fig):
    os.makedirs(path_fig)
filobs = path_out+'observation_MIR' # after input_presimulation.py

## Read h5 file
spec_unit = read_hdf5(filobs, 'spectral unit')[0]
wvl = read_hdf5(filobs, 'wavelength (microns)')
FnuOBS = read_hdf5(filobs, 'FnuOBS ('+spec_unit+')')
dFnuOBS = read_hdf5(filobs, 'dFnuOBS ('+spec_unit+')')
FnuMOD = read_hdf5(filout, 'FnuMOD ('+spec_unit+')')
FnuCONT = read_hdf5(filout, 'FnuCONT ('+spec_unit+')')
FnuLINE = read_hdf5(filout, 'FnuLINE ('+spec_unit+')')
FnuBAND = read_hdf5(filout, 'FnuBAND ('+spec_unit+')')
FnuSTAR = read_hdf5(filout, 'FnuSTAR ('+spec_unit+')')
Pabs = read_hdf5(filout, 'PABS')

## Make plot tables
Ncont = FnuCONT.shape[0]
Nline = FnuLINE.shape[0]
Nband = FnuBAND.shape[0]
Nstar = FnuSTAR.shape[0]
Fnu_tab = []
lab_tab = []
col_tab = ['k'] # pplot() none data iter
Fnu_tab.append(FnuOBS)
lab_tab.append('Obs')
col_tab.append('y')
Fnu_tab.append(FnuMOD)
lab_tab.append('Total')
col_tab.append('k')
for i in range(Ncont):
    Fnu_tab.append(FnuCONT[i,:,:,:]*np.prod(Pabs,axis=0))
    if i==0:
        lab_tab.append('cont')
    else:
    	lab_tab.append('')
    col_tab.append('pink')
for i in range(Nline):
    Fnu_tab.append((FnuLINE[i,:,:,:]+np.sum(FnuCONT,axis=0)+FnuSTAR[0])*np.prod(Pabs,axis=0))
    if i==0:
        lab_tab.append('line')
    else:
    	lab_tab.append('')
    col_tab.append('g')
for i in range(Nband):
    Fnu_tab.append((FnuBAND[i,:,:,:]+np.sum(FnuCONT,axis=0)+FnuSTAR[0])*np.prod(Pabs,axis=0))
    if i==0:
        lab_tab.append('band')
    else:
    	lab_tab.append('')
    col_tab.append('b')
for i in range(Nstar):
    Fnu_tab.append(FnuSTAR[i,:,:,:]*np.prod(Pabs,axis=0))
    if i==0:
        lab_tab.append('star')
    else:
    	lab_tab.append('')
    col_tab.append('orange')

## Plot individual fit
##---------------------
filename = path_fig+'presimulation'

xlab = r'$\lambda\,(\mu m)$'
if spec_unit=='MKS':
    ylab = r'Surface brightness $\,(W/m2/Hz/sr)$'
elif spec_unit=='MJyovsr':
    ylab=r'Surface brightness $\,(MJy/sr)$'
Fnu_tab = np.array(Fnu_tab)

p = pplot(xlog=1, ylog=1, nonposx='clip', nonposy='clip',
          xlim=(2.5, 22.),
          ylim=(1e0, 1e3),
          xlabel=xlab, ylabel=ylab,
          legend='upper left', clib=col_tab)
for i in range(Fnu_tab.shape[0]):
    p.add_plot(wvl, Fnu_tab[i,:,0,0], label=lab_tab[i])
    
# p = pplot(wvl, FnuOBS/dFnuOBS, 
#           xlog=0, ylog=0, xlim=(4.5, 22.), #ylim=(0,2.), 
#           xlab=xlab, ylab=ylab+'/SN ratio', 
#           legend='best', label='S/N', c='r')
# p.add_plot(wvl, FnuOBS, label='Data', c='k')

p.save(filename)
# plt.show()

print('>>> Coucou show_presimulation= [Done] <<<')

