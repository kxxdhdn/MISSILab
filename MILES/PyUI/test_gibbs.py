#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""

This is the visualization of the test_gibbs

"""

import os
import math
import numpy as np
import matplotlib.pyplot as plt
# from astropy.constants import M_sun, pc, 
# from astropy.units import Jy,

## astylo
from astylo.iolib import read_hdf5

## Path
mroot = os.path.dirname(os.path.abspath(__file__))+'/../'
path_out = mroot+'tests/out/'
filename = 'test_gibbs'
flogname = 'parlog_gibbs'

## Read h5 files
wvl, FnuOBS, FnuMOD = \
    read_hdf5(path_out+filename, 
    'Wavelength (microns)', 'FnuOBS (MJyovsr)', 'FnuMOD (MJyovsr)')
Fnu1par, lnLHobs, testgrid = read_hdf5(path_out+filename, 
    'var1par', 'lnLHobs', 'testgrid')
parmu, parsig = read_hdf5(path_out+filename, 
    'Mean of parameter value', 'Sigma of parameter value')

Nmcmc, pararr = read_hdf5(path_out+flogname, 
	'Length of MCMC', 'Parameter values')

## Vary one param
##----------------
# y = [FnuOBS]
# y.extend(Fnu1par)
# p1 = plot2d_m(wvl, y, 
#              xall=wvl, xlog=1, ylog=1, xlim=(1.,40.), ylim=(0.,2.), 
#              xlab=r'$\lambda\,(micron)$', ylab=r'$F_{\nu} \,(MJy/sr)$', 
#              lablist=[str(i) for i in range(6)], 
#              legend='best', cl=['y', 'r', 'orange', 'g', 'b', 'purple'])
# p1.save(path_out+filename+'_var1par.png')

## Plot
##======
ipar = 1
Nbin = int(Nmcmc)
fig, axes = plt.subplots(nrows=2, ncols=1, figsize=(8,9))
plt.subplots_adjust(left=.1, bottom=.05, \
                    right=.99, top=.95, wspace=0., hspace=.2)
ax1, ax2 = axes
## Param dist
##------------
n, bins, patches = ax1.hist(pararr[:,ipar-1], Nbin, density=1, facecolor='g', alpha=.5)
# print(pararr[:,ipar-1])
ax1.set_xlabel('param value')
ax1.set_ylabel('pdf')
ax1.set_title(fr'$Param\ dist\ \mu={parmu[ipar-1]},\ \sigma={parsig[ipar-1]}$')
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
# ax2.set_xlim((3.,40.))
ax2.set_ylim((1.E-2,2.))
ax2.set_xlabel(r'$\lambda\,(micron)$')
ax2.set_ylabel(r'$F_{\nu} \,(MJy/sr)$')
ax2.legend()
ax2.set_title('Spectral fitting')

fig.savefig(path_out+filename+'.png')
plt.show()

print('>>> Coucou test_gibbs [Done] <<<')
