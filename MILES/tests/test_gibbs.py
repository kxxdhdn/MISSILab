#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""

This is the visualization of test_gibbs

"""

import os
import math
import numpy as np
import matplotlib.pyplot as plt
# from astropy.constants import M_sun, pc, 
# from astropy.units import Jy,

## rapyuta
from rapyuta.inout import read_hdf5

## Path
mroot = os.path.dirname(os.path.abspath(__file__))+'/'
path_out = mroot+'out/'
filename = 'test_gibbs'
flogname = 'parlog_gibbs'

## Read h5 files
wvl = read_hdf5(path_out+filename, 'Wavelength (microns)')
FnuOBS = read_hdf5(path_out+filename, 'FnuOBS (MJyovsr)')
FnuMOD = read_hdf5(path_out+filename, 'FnuMOD (MJyovsr)')
pargen = read_hdf5(path_out+filename, 'True parameter value')
# Fnu1par = read_hdf5(path_out+filename, 'var1par')
# lnLHobs = read_hdf5(path_out+filename, 'lnLHobs')
# testgrid = read_hdf5(path_out+filename, 'testgrid')
parmu = read_hdf5(path_out+filename, 'Mean of parameter value')
parsig = read_hdf5(path_out+filename, 'Sigma of parameter value')
t_burnin = read_hdf5(path_out+filename, 't_burnin')[0]
t_end = read_hdf5(path_out+filename, 't_end')[0]

Nmcmc = read_hdf5(path_out+flogname, 'Length of MCMC')
pararr = read_hdf5(path_out+flogname, 'Parameter values')

## Plot
##======
ipar = 7
Nbin = int(Nmcmc)
fig, axes = plt.subplots(nrows=2, ncols=1, figsize=(8,8))
plt.subplots_adjust(left=.1, bottom=.1, \
                    right=.95, top=.95, wspace=0., hspace=.3)
ax1, ax2 = axes
## Param dist
##------------
n, bins, patches = ax1.hist(pararr[t_burnin:t_end,ipar], Nbin, density=1, facecolor='g', alpha=.5)
# print(pararr[:,ipar])
ax1.axvline(x=pargen[ipar], color='r', lw=2., label='True')
ax1.axvline(x=parmu[ipar], color='b', label='Mean')
ax1.axvline(x=parmu[ipar]-3*parsig[ipar], color='b', ls='--', label=r'$-3\sigma$')
ax1.axvline(x=parmu[ipar]+3*parsig[ipar], color='b', ls='--', label=r'$+3\sigma$')
# ax1.set_xscale('log')
ax1.set_xlabel('param value')
ax1.set_ylabel('pdf')
ax1.legend()
ax1.set_title(fr'$Parameter\ (true: {pargen[ipar]:.2f})\ distribution\ \mu={parmu[ipar]:.2f},\ \sigma={parsig[ipar]:.2f}$')
## LH dist
##---------
# ax2.plot(testgrid, lnLHobs)
# ax2.set_xlabel('pargrid')
# ax2.set_ylabel('lnLHobs')
# ax2.set_title('lnLHobs')
## Mod vs Obs
##------------
ax2.plot(wvl, FnuOBS, c='y', label='Obs')
ax2.plot(wvl, FnuMOD, c='k', label='Mod')
ax2.set_xscale('log')
ax2.set_yscale('log')
ax2.set_xlim((2.,41.))
# ax2.set_ylim((1.E-2,2.))
ax2.set_xlabel(r'$\lambda\,(\mu m)$')
ax2.set_ylabel(r'$F_{\nu} \,(MJy/sr)$')
# ax2.set_ylabel(r'$F_{\nu} \,(W/m2/Hz/sr)$')
ax2.legend()
ax2.set_title('Spectral fitting')

fig.savefig(path_out+filename+'.png')
# plt.show()

print('>>> Coucou test_gibbs [Done] <<<')
