#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""

This is the visualization of test_fitHB

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
filename = 'test_fitHB'
flogname = 'parlog_fitHB'

x = 1
y = 1
## Read h5 file
spec_unit = read_hdf5(filobs, 'Spectral unit')[0]
wvl = read_hdf5(filobs, 'Wavelength (microns)')
FnuOBS = read_hdf5(filobs, 'FnuOBS ('+spec_unit+')')[:,y,x]
dFnuOBS = read_hdf5(filobs, 'dFnuOBS ('+spec_unit+')')[:,y,x]
FnuMOD = read_hdf5(path_out+filename, 'FnuMOD ('+spec_unit+')')[:,y,x]
FnuCONT = read_hdf5(path_out+filename, 'FnuCONT ('+spec_unit+')')[:,:,y,x]
FnuLINE = read_hdf5(path_out+filename, 'FnuLINE ('+spec_unit+')')[:,:,y,x]
FnuBAND = read_hdf5(path_out+filename, 'FnuBAND ('+spec_unit+')')[:,:,y,x]
FnuSTAR = read_hdf5(path_out+filename, 'FnuSTAR ('+spec_unit+')')[:,:,y,x]
parmu = read_hdf5(path_out+filename, 'Mean of parameter value')[:,y,x]
parsig = read_hdf5(path_out+filename, 'Sigma of parameter value')[:,y,x]

Nmcmc = read_hdf5(path_out+flogname, 'Length of MCMC')
parnam = read_hdf5(path_out+flogname, 'Parameter label')
pararr = read_hdf5(path_out+flogname, 'Parameter values')[:,:,y,x]
print('Nmcmc', Nmcmc)
# print(pararr.shape)
# print(FnuCONT)
# print(FnuLINE)
# print(FnuBAND)
# print(FnuSTAR)
# exit()

## Make plot tables
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
    Fnu_tab.append(FnuCONT[i,:])
    if i==0:
        lab_tab.append('cont')
    else:
    	lab_tab.append('')
    col_tab.append('pink')
for i in range(Nline):
    Fnu_tab.append(FnuLINE[i,:])
    if i==0:
        lab_tab.append('line')
    else:
    	lab_tab.append('')
    col_tab.append('g')
for i in range(Nband):
    Fnu_tab.append(FnuBAND[i,:])
    if i==0:
        lab_tab.append('band')
    else:
    	lab_tab.append('')
    col_tab.append('b')
for i in range(Nstar):
    Fnu_tab.append(FnuSTAR[i,:])
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

## Plot
##======
ipar = 3
Nbin = int(Nmcmc)
fig, axes = plt.subplots(nrows=2, ncols=1, figsize=(8,9))
plt.subplots_adjust(left=.1, bottom=.05, \
                    right=.99, top=.95, wspace=0., hspace=.2)
ax1, ax2 = axes
## Param dist
##------------
# print(pararr[:,ipar-1].shape)
n, bins, patches = ax1.hist(pararr[:,ipar-1], Nbin, density=1, facecolor='g', alpha=.5)
ax1.set_xlabel('param value')
ax1.set_ylabel('pdf')
ax1.set_title(fr'${parnam[ipar-1]}\ dist\ \mu={parmu[ipar-1]},\ \sigma={parsig[ipar-1]}$')
## LH dist
##---------
# ax2.plot(testgrid, lnLHobs)
# ax2.set_xlabel('pargrid')
# ax2.set_ylabel('lnLHobs')
# ax2.set_title('lnLHobs')

## Show fit
##----------
p = plot2d_m(wvl, Fnu_tab, 
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

print('>>> Coucou test_fitHB [Done] <<<')

