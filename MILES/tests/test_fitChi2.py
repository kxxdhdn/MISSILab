#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""

This is the visualization of test_fitChi2

"""
import os
import numpy as np
import matplotlib.pyplot as plt

## astylo
from astylo.iolib import read_hdf5
from astylo.plib import plot2d_m

## Path
mroot = os.path.dirname(os.path.abspath(__file__))+'/'
filobs = mroot+'dat/observations_fitMIR'
path_out = mroot+'out/'
filename = 'test_fitChi2'

# x = 46
# y = 20
## Read h5 file
spec_unit = read_hdf5(filobs, 'Spectral unit')[0]
wvl = read_hdf5(filobs, 'Wavelength (microns)')
FnuOBS = read_hdf5(filobs, 'FnuOBS ('+spec_unit+')')#[:,y,x]
dFnuOBS = read_hdf5(filobs, 'dFnuOBS ('+spec_unit+')')#[:,y,x]
FnuMOD = read_hdf5(path_out+filename, 'FnuMOD ('+spec_unit+')')#[:,y,x]
FnuCONT = read_hdf5(path_out+filename, 'FnuCONT ('+spec_unit+')')#[:,:,y,x]
FnuLINE = read_hdf5(path_out+filename, 'FnuLINE ('+spec_unit+')')#[:,:,y,x]
FnuBAND = read_hdf5(path_out+filename, 'FnuBAND ('+spec_unit+')')#[:,:,y,x]
FnuSTAR = read_hdf5(path_out+filename, 'FnuSTAR ('+spec_unit+')')#[:,:,y,x]

## Make plot tables
Nx = FnuOBS.shape[2]
Ny = FnuOBS.shape[1]
Ncont = FnuCONT.shape[0]
Nline = FnuLINE.shape[0]
Nband = FnuBAND.shape[0]
Nstar = FnuSTAR.shape[0]
Fnu_tab = []
lab_tab = []
col_tab = []
Fnu_tab.append(FnuOBS)
lab_tab.append('Obs')
col_tab.append('y')
Fnu_tab.append(FnuMOD)
lab_tab.append('Total')
col_tab.append('k')
for i in range(Ncont):
    Fnu_tab.append(FnuCONT[i,:,:,:])
    if i==0:
        lab_tab.append('cont')
    else:
    	lab_tab.append('')
    col_tab.append('pink')
for i in range(Nline):
    Fnu_tab.append(FnuLINE[i,:,:,:])
    if i==0:
        lab_tab.append('line')
    else:
    	lab_tab.append('')
    col_tab.append('g')
for i in range(Nband):
    Fnu_tab.append(FnuBAND[i,:,:,:])
    if i==0:
        lab_tab.append('band')
    else:
    	lab_tab.append('')
    col_tab.append('b')
for i in range(Nstar):
    Fnu_tab.append(FnuSTAR[i,:,:,:])
    if i==0:
        lab_tab.append('star')
    else:
    	lab_tab.append('')
    col_tab.append('orange')

xlab = r'$\lambda\,(\mu m)$'
if spec_unit=='MKS':
    ylab = r'Surface brightness $\,(W/m2/Hz/sr)$'
elif spec_unit=='MJyovsr':
    ylab=r'Surface brightness $\,(MJy/sr)$'
Fnu_tab = np.array(Fnu_tab)

# for x in range(Nx):
#     for y in range(Ny):
p = plot2d_m(wvl, Fnu_tab[:,:,1,0], 
             xall=wvl, xlog=0, ylog=0, xlim=(4.5, 22.), #ylim=(0.,1.e3), 
             xlab=xlab, ylab=ylab, 
             lablist=lab_tab, 
             legend='best', cl=col_tab)
# p = plot2d_m(wvl, [FnuOBS/dFnuOBS, FnuOBS], 
#              xall=wvl, xlog=0, ylog=0, xlim=(4.5, 22.), #ylim=(0,2.), 
#              xlab=xlab, ylab=ylab+'/SN ratio', 
#              lablist=['S/N', 'Data'], 
#              legend='best', cl=['r', 'k'])
p.save(path_out+filename+'.png')

plt.show()

print('>>> Coucou test_fitChi2 [Done] <<<')

