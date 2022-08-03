#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""

This is the visualization of simulation_MIR.h5 (simulated spectra)

"""

import os
import numpy as np
from matplotlib.ticker import ScalarFormatter, NullFormatter

## rapyuta
from rapyuta.arrays import closest
from rapyuta.inout import read_hdf5
from rapyuta.plots import pplot

## local
from auxil import croot

clist = ['k',
         'c',
         'g',
         'm',
         'orange']

labelist = ['c04_SN2',
            'c20_SN2',
            'c04_SN100',
            'c20_SN100',]

## Path
##------
path_out = croot+'/../out/'
path_fig = path_out+'Figures/'
if not os.path.exists(path_fig):
    os.makedirs(path_fig)
filout = path_out+'simulation_MIR'

## Read h5 file
spec_unit = read_hdf5(filout, 'spectral unit')[0]
wvl = read_hdf5(filout, 'wavelength (microns)')
FnuOBS = read_hdf5(filout, 'FnuOBS ('+spec_unit+')')
dFnuOBS = read_hdf5(filout, 'dFnuOBS ('+spec_unit+')')

## Added a NaN at 5 micron
i5 = closest(wvl, 5., 'left') + 1
wvl = np.insert(wvl, i5, 5.)
FnuOBS = np.insert(FnuOBS, i5, np.nan, axis=0)
dFnuOBS = np.insert(dFnuOBS, i5, np.nan, axis=0)

Nw, Ny, Nx = FnuOBS.shape

## Negative to NaN
FnuOBS[FnuOBS<0] = np.nan

## Make plots
xlab = r'$\lambda\,(\mu m)$'
if spec_unit=='MKS':
    ylab = r'$F_{\nu}\ \,(W/m2/Hz/sr)$'
elif spec_unit=='MJyovsr':
    ylab = r'$F_{\nu}\ \,(MJy/sr)$'

## Create data tables
lab_tab = np.empty((Ny,Nx), dtype=('<U30'))
titlist = np.empty((Ny,Nx), dtype=('<U60'))
Fnu_tab = []
dFnu_tab = []
for x in range(Nx):
    for y in range(Ny):
        lab_tab[y,x] = labelist[y]+'_'+str(x+1)
        titlist[y,x] = 'Synthetic Spectra '+lab_tab[y,x]+' (rest frame)'
        Fnu_tab.append(FnuOBS[:,y,x])
        dFnu_tab.append(dFnuOBS[:,y,x])
        
Fnu_tab = np.array(Fnu_tab).reshape((Nx,Ny,Nw))
dFnu_tab = np.array(dFnu_tab).reshape((Nx,Ny,Nw))

## Plot individual spectrum with S/N
##------------------------------------
for y in range(Ny):
    for x in range(Nx):
        filename1 = path_fig+'simu_'+str(x+1)+'_'+str(y+1)+'.png'

        p = pplot(xlog=1, ylog=1,# nonposx='clip', nonposy='clip',
                  xlim=(2.4,21), ylim=(1e-1,1e4),
                  xlabel=xlab, ylabel=ylab,
                  title=titlist[y,x],
                  figsize=(10,8), left=.15, right=.95, top=.9, bottom=.1,
                  titlesize=20, labelsize=20, ticksize=20)
        p.add_plot(wvl, Fnu_tab[x,y,:], yerr=dFnu_tab[x,y,:],
                   ec='grey', elw=1, c='k', label='Data')
        p.add_plot(wvl, dFnu_tab[x,y,:], c='r', label='Errors')

        # p.ax.text(.9,.05,'(a)',size=20,c='grey',transform=p.ax.transAxes)
        xtic = [2.5, 3, 3.5, 4, 5, 6, 7, 8, 10, 12, 15, 20,]# 30, 40]
        xtic_min = np.arange(2.5, 20., .1)
        p.ax.set_xticks(xtic, minor=False) # major
        p.ax.set_xticks(xtic_min, minor=True) # minor
        p.ax.xaxis.set_major_formatter(ScalarFormatter()) # major
        p.ax.xaxis.set_minor_formatter(NullFormatter()) # minor
        p.ax.legend(loc='upper left',# bbox_to_anchor=(1,1),
                    fontsize=20, framealpha=0)
                    
        p.save(filename1, transparent=True)

## Plot stacked spectra
##----------------------
for x in range(Nx):
    filename2 = path_fig+'simu_stack_x='+str(x+1)+'.png'

    p = pplot(xlog=1, ylog=1, nonposx='clip', nonposy='clip',
              xlim=(2.4,21), ylim=(1e0,2e4),
              xlabel=xlab, ylabel=ylab,
              figsize=(10,8), left=.15, right=.95, top=.95, bottom=.1,
              titlesize=20, labelsize=20, ticksize=20, title=None, clib=clist)
    for y in range(Ny):
        p.add_plot(wvl, Fnu_tab[x,y,:], label=lab_tab[y,x])

    # p.ax.text(.9,.05,'(a)',size=20,c='grey',transform=p.ax.transAxes)
    xtic = [2.5, 3, 3.5, 4, 5, 6, 7, 8, 10, 12, 15, 20,]# 30, 40]
    xtic_min = np.arange(2.5, 20., .1)
    p.ax.set_xticks(xtic, minor=False) # major
    p.ax.set_xticks(xtic_min, minor=True) # minor
    p.ax.xaxis.set_major_formatter(ScalarFormatter()) # major
    p.ax.xaxis.set_minor_formatter(NullFormatter()) # minor
    p.ax.legend(loc='upper left',# bbox_to_anchor=(1,1),
                fontsize=20, framealpha=0)

    p.save(filename2, transparent=True)
    
                
print('>>> Coucou show_simulation [Done] <<<')
