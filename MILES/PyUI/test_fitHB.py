#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""

This is the visualization of test_fitHBsyn

"""

import os
import math
import numpy as np
import matplotlib.pyplot as plt

## astylo
from astylo.iolib import read_hdf5
from astylo.plib import plot2d_m, plot2d

## Path
mroot = os.path.dirname(os.path.abspath(__file__))+'/../'
filobs = mroot+'tests/dat/observations_fitMIR'
path_out = mroot+'tests/out/'
filename = 'test_fitHB'
flogname = 'parlog_fitHB'

## Read h5 file
wvl = read_hdf5(filobs, 'Wavelength (microns)')
FnuOBS = read_hdf5(filobs, 'FnuOBS (MKS)')[:,6,5]
FnuMOD = read_hdf5(path_out+filename, 'FnuMOD (MKS)')[:,6,5]
FnuCONT = read_hdf5(path_out+filename, 'FnuCONT (MKS)')[:,6,5]
FnuLINE = read_hdf5(path_out+filename, 'FnuLINE (MKS)')[:,6,5]
FnuBAND = read_hdf5(path_out+filename, 'FnuBAND (MKS)')[:,6,5]
FnuSTAR = read_hdf5(path_out+filename, 'FnuSTAR (MKS)')[:,6,5]
parmu = read_hdf5(path_out+filename, 'Mean of parameter value')[:,6,5]
parsig = read_hdf5(path_out+filename, 'Sigma of parameter value')[:,6,5]

Nmcmc = read_hdf5(path_out+flogname, 'Length of MCMC')
pararr = read_hdf5(path_out+flogname, 'Parameter values')

## Param & LH dist
##-----------------
ipar = 7
Nbin = int(Nmcmc/10)
fig, axes = plt.subplots(nrows=2, ncols=1, figsize=(8,9))
plt.subplots_adjust(left=.1, bottom=.05, \
                    right=.99, top=.95, wspace=0., hspace=.2)
ax1, ax2 = axes
## Param
n, bins, patches = ax1.hist(pararr[ipar-1], Nbin, density=1, facecolor='g', alpha=.5)
# print(pararr[:,ipar-1])
ax1.set_xlabel('param value')
ax1.set_ylabel('pdf')
ax1.set_title(fr'$Param\ dist\ \mu={parmu[ipar-1]},\ \sigma={parsig[ipar-1]}$')
## Spectral fitting
##------------------
ax2.plot(wvl, FnuOBS, c='y', label='Obs')
ax2.plot(wvl, FnuMOD, c='k', label='Total')
ax2.plot(wvl, FnuCONT, c='pink', label='cont')
ax2.plot(wvl, FnuLINE, c='g', label='line')
ax2.plot(wvl, FnuBAND, c='b', label='band')
ax2.plot(wvl, FnuSTAR, c='orange', label='star')
ax2.set_xscale('log')
ax2.set_yscale('log')
# ax2.set_xlim((2.,40.))
ax2.set_ylim((1.E-2,2.))
ax2.set_xlabel(r'$\lambda\,(micron)$')
ax2.set_ylabel(r'$F_{\nu} \,(MJy/sr)$')
ax2.legend('best')
ax2.set_title('Spectral fitting')

# p = plot2d_m(wvl, [FnuOBS, FnuMOD, FnuCONT, FnuLINE, FnuBAND, FnuSTAR], 
#              xall=wvl, xlog=1, ylog=1, xlim=(2., 40.), ylim=(0,2.), 
#              xlab=r'$\lambda\,(micron)$', ylab=r'$F_{\nu} \,(MJy/sr)$', 
#              lablist=['Obs', 'Total', 'cont', 'lines', 'bands', 'star'], 
#              legend='best', cl=['y', 'k', 'pink', 'g', 'b', 'orange'])
# p.save(path_out+filename+'.png')

fig.savefig(path_out+filename+'.png')
plt.show()

print('>>> Coucou test_fitHBsyn [Done] <<<')
