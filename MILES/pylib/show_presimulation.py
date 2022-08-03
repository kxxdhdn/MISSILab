#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""

This is the visualization of presimulation.h5 (chi2 fit)

"""
import os
import numpy as np
import matplotlib.pyplot as plt

## laputan
from laputan.arrays import closest
from laputan.inout import read_hdf5
from laputan.plots import pplot

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
par = read_hdf5(filout, 'Best fitted parameter value')
parname = read_hdf5(filout, 'Best fitted parameter label')

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

p = pplot(xlog=0, ylog=0, xlim=(4., 22.), #ylim=(0.,1.e3), 
          xlab=xlab, ylab=ylab, 
          legend='upper left', clib=col_tab)
for i in range(Fnu_tab.shape[0]):
    p.add_plot(wvl, Fnu_tab[i,:,0,0], label=lab_tab[i])
    
# p = pplot(wvl, FnuOBS/dFnuOBS, 
#           xlog=0, ylog=0, xlim=(4.5, 22.), #ylim=(0,2.), 
#           xlab=xlab, ylab=ylab+'/SN ratio', 
#           legend='best', label='S/N', c='r')
# p.add_plot(wvl, FnuOBS, label='Data', c='k')

p.save(filename)

## Plot individual band in table
##-------------------------------
# Cband = []
# for b in range(Nband):
#     ipar = Ncont*2+Nline*3+b*4+1
#     Cband.append(par[ipar])
#     # print(parname[ipar])
# Cband = np.array(Cband)

# for k in range(int(Nband/9)+1):
#     filename1 = path_fig+'presimulation_band'+str(k*9+1)+'-'+str(k*9+9)+'.png'

#     fig, axes = plt.subplots(nrows=3, ncols=3, figsize=(9,9))
#     plt.subplots_adjust(left=.1, bottom=.05, \
#                         right=.99, top=.95, wspace=.3, hspace=.4)
#     for j in range(9):
#         b = 9*k+j
#         if b<Nband:
#             valband = FnuBAND[b,closest(wvl,Cband[b]),0,0]
#             px, py = int(j/3), j%3
#             # axes[px,py].plot(wvl, Fnu_tab[0,:,0,0],c='k')
#             axes[px,py].plot(wvl, FnuBAND[b,:,0,0],c='b')
#             axes[px,py].vlines(Cband[b],0,1.5*valband,colors='r')
#             axes[px,py].set_xlim((Cband[b]-.2, Cband[b]+.2))
#             axes[px,py].set_ylim((0, 1.5*valband))
#             axes[px,py].set_xlabel(xlab)
#             axes[px,py].set_ylabel(ylab)
#             axes[px,py].set_title(parname[Ncont*2+Nline*3+b*4+1])
    
#     plt.savefig(filename1)

# plt.show()

print('>>> Coucou show_presimulation= [Done] <<<')

