#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""

MIR inter-calibration scheme

"""

# import warnings
# warnings.filterwarnings("ignore", category=RuntimeWarning)


import numpy as np
from matplotlib.ticker import ScalarFormatter, NullFormatter

## rapyuta
from rapyuta.inout import fclean, fitsext, read_fits, write_fits, read_hdf5
from rapyuta.arrays import closest
from rapyuta.plots import pplot

## Local
from buildinfo import ( src, path_cal, path_fig, path_out, parobs, )


## Banner
print('\n============================================================\n')

print('        MIRAGE - MIR spectrum inter-calibration - '+src)

print('\n============================================================\n')


x, y = 0, 7 # x=0, y=0-41
iobs = 0 # 0-13
Name = ['A','B','C','D','E','F','G','H','I','J','K','L','M','N']
ds = read_fits(path_out+src+'_'+parobs[iobs][0], path_out+src+'_'+parobs[iobs][0]+'_unc')
ds1 = read_fits(path_out+src+'_'+parobs[iobs][0]+'_IRC', path_out+src+'_'+parobs[iobs][0]+'_IRC_unc')
ds2 = read_fits(path_out+src+'_'+parobs[iobs][0]+'_IRS', path_out+src+'_'+parobs[iobs][0]+'_IRS_unc')
irac1 = read_fits(path_cal+src+'_IRAC1_SINGS_'+parobs[iobs][0], path_cal+src+'_IRAC1_SINGS_'+parobs[iobs][0]+'_unc')
irc = read_fits(path_cal+src+'_IRAC1_IRC_'+parobs[iobs][0], path_cal+src+'_IRAC1_IRC_'+parobs[iobs][0]+'_unc')
irac4 = read_fits(path_cal+src+'_IRAC4_SINGS_'+parobs[iobs][0], path_cal+src+'_IRAC4_SINGS_'+parobs[iobs][0]+'_unc')
sl = read_fits(path_cal+src+'_IRAC4_IRS_'+parobs[iobs][0], path_cal+src+'_IRAC4_IRS_'+parobs[iobs][0]+'_unc')
mips1 = read_fits(path_cal+src+'_MIPS1_SINGS_'+parobs[iobs][0], path_cal+src+'_MIPS1_SINGS_'+parobs[iobs][0]+'_unc')
ll = read_fits(path_cal+src+'_MIPS1_IRS_'+parobs[iobs][0], path_cal+src+'_MIPS1_IRS_'+parobs[iobs][0]+'_unc')
# Nw, Ny, Nx = ds.data.shape
# for x in range(1):
#     for y in range(4,5):
#         for k in range (Nw):
i0 = closest(ds.wave, 5) + 1
i1 = closest(ds.wave, 14.29, 'left') + 1
g0 = irac1.data[y,x]/irc.data[y,x]
g1 = irac4.data[y,x]/sl.data[y,x]
g2 = mips1.data[y,x]/ll.data[y,x]
ds.data[:i0,y,x] = ds.data[:i0,y,x] * g0
ds.data[i0:i1,y,x] = ds.data[i0:i1,y,x] * 1.05#g1
ds.data[i1:,y,x] = ds.data[i1:,y,x] * g2
print(g0, g1, g2)
p = pplot(ds.wave, ds.data[:,y,x], yerr=ds.unc[:,y,x],
          xlog=1, ylog=1, 
          lw=1, c='k', ec='r', label='Combined calib',
          xlim=(1.75,40),
          # ylim=ylim,
          xlabel=r'${\rm Wavelengths,}\ \lambda\ (\mu m)$',
          ylabel=r'$\rm F_{\nu}\ (MJy/sr)$',
          figsize=(16,8),# top=.99, right=.85, left=.1,
          loc='upper left', legendalpha=0,
          xysize=15, tksize=15, legendsize=15)
p.add_plot(ds1.wave, ds1.data[:,y,x], ls='dashed', c='grey', label='IRC no calib')
p.add_plot(ds2.wave, ds2.data[:,y,x], ls='dashed', c='grey', label='IRS no calib')
p.add_plot(3.6, irac1.data[y,x], marker='o', ms=15, c='r', label=r'$\rm IRAC_{3.6\mu m}$')
p.add_plot(3.6, irc.data[y,x], marker='^', ms=15, c='r', alpha=0.5, label=r'$\rm sIRAC_{3.6\mu m}$')
p.add_plot(8, irac4.data[y,x], marker='o', ms=15, c='g', label=r'$\rm IRAC_{8.0\mu m}$')
p.add_plot(8, irac4.data[y,x]/1.05, marker='^', ms=15, c='g', alpha=0.5, label=r'$\rm sIRAC_{8.0\mu m}$')
p.add_plot(24, mips1.data[y,x], marker='o', ms=15, c='c', label=r'$\rm MIPS_{24\mu m}$')
p.add_plot(24, ll.data[y,x], marker='^', ms=15, c='c', alpha=0.5, label=r'$\rm sMIPS_{24\mu m}$')
xtic = [2, 3, 4, 5, 6, 7, 8, 10, 12, 15, 20, 30, 40]
xtic_min = np.arange(2., 40., .5)
p.ax.set_xticks(xtic, minor=False) # major
p.ax.set_xticks(xtic_min, minor=True) # minor
p.ax.xaxis.set_major_formatter(ScalarFormatter()) # major
p.ax.xaxis.set_minor_formatter(NullFormatter()) # minor
p.ax.xaxis.set_ticks_position('both') # ticks
p.ax.yaxis.set_ticks_position('both')
p.ax.tick_params(axis='both', which='both', direction='in')
p.ax.tick_params(which='major', length=8, width=2)
p.ax.tick_params(which='minor', length=4, width=2)

p.save(path_cal+src+'_IC_'+Name[iobs]+str(int(y/parobs[iobs][6])), figtight=True)

