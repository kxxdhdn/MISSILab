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

labelist = ['c1 SN10',
            'c10 SN10',
            'c1 SN100',
            'c10 SN100',]

## Path
##------
path_out = croot+'/../out/'
path_fig = path_out+'Figures/'
os.makedirs(path_fig, exist_ok=True)
filout = path_out+'simulation_MIR'

## Read h5 file
spec_unit = read_hdf5(filout, 'spectral unit')[0]
wvl = read_hdf5(filout, 'wavelength (microns)')
FnuOBS = read_hdf5(filout, 'FnuOBS ('+spec_unit+')') # FnuSIM
dFnuOBS = read_hdf5(filout, 'dFnuOBS ('+spec_unit+')') # dFnuSIM

## Added a NaN at 5 micron
i5 = closest(wvl, 5., 'left')
FnuOBS[i5,:,:] = np.nan
dFnuOBS[i5,:,:] = np.nan
# wvl = np.insert(wvl, i5, 5.)
# FnuOBS = np.insert(FnuOBS, i5, np.nan, axis=0)
# dFnuOBS = np.insert(dFnuOBS, i5, np.nan, axis=0)

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
        lab_tab[y,x] = labelist[y]+' #'+str(x+1)
        titlist[y,x] = 'Synthetic Spectra '+lab_tab[y,x]+' (rest frame)'
        Fnu_tab.append(FnuOBS[:,y,x])
        dFnu_tab.append(dFnuOBS[:,y,x])
        
Fnu_tab = np.array(Fnu_tab).reshape((Nx,Ny,Nw))
dFnu_tab = np.array(dFnu_tab).reshape((Nx,Ny,Nw))

xtic = [2.5, 3, 4, 5, 7, 9, 12, 15, 20,]# 30, 40]
xtic_min = np.arange(2.5, 20., .1)

## Plot individual spectrum with S/N
##------------------------------------
for y in range(Ny):
    for x in range(Nx):

        p = pplot(xlog=1, ylog=1,# nonposx='clip', nonposy='clip',
                  ylim=(1e-1,5e3),
                  xlabel=xlab, ylabel=ylab,
                  title=titlist[y,x],
                  xtk=xtic, xtkmi=xtic_min, xtkform='mylog', ytkform='log_sci',
                  figsize=(10,8), loc='upper left', legendalpha=0,
                  titlesize=20, xysize=20, tksize=20, legendsize=20)
        p.add_plot(wvl, Fnu_tab[x,y,:], yerr=dFnu_tab[x,y,:],
                   ec='grey', elw=1, c='k', label='Data')
        p.add_plot(wvl, dFnu_tab[x,y,:], c='r', label='Errors')
        # p.ax.text(.9,.05,'(a)',size=20,c='grey',transform=p.ax.transAxes)
                    
        filename = path_fig+'simu_'+str(x+1)+'_'+str(y+1)+'.png'
        p.save(filename, transparent=True, figtight=True)

## Plot stacked spectra
##----------------------
for x in range(Nx):
    p = pplot(xlog=1, ylog=1,# nonposx='clip', nonposy='clip',
              ylim=(2e-1,2e4),
              xlabel=xlab, ylabel=ylab, clib=clist,
              xtk=xtic, xtkmi=xtic_min, xtkform='mylog', ytkform='log_sci',
              figsize=(10,8), loc='upper left', legendalpha=0,
              titlesize=20, xysize=20, tksize=20, legendsize=20)
    for y in range(Ny):
        ## scale by 0.1 to separate spectra]
        if y==0 or y==2:
            Fnu_tab[x,y,:] *= 0.1
            lab_tab[y,x] += ' (scaled by 0.1)'
        p.add_plot(wvl, Fnu_tab[x,y,:], label=lab_tab[y,x])
    # p.ax.text(.9,.05,'(a)',size=20,c='grey',transform=p.ax.transAxes)

    filename = path_fig+'simu_stack_x='+str(x+1)+'.png'
    p.save(filename, transparent=True, figtight=True)
    
                
print('>>> Coucou show_simulation [Done] <<<')
