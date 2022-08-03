#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""

This is the visualization of presimulation.h5 (chi2 fit)

"""

import os, pathlib
import numpy as np
from matplotlib.ticker import ScalarFormatter, NullFormatter

## rapyuta
from rapyuta.arrays import closest
from rapyuta.inout import read_hdf5
from rapyuta.plots import pplot

## local
from auxil import croot, calexpansion
                       
## Path
##------
path_out = croot+'/../out/'
path_fig = path_out+'Figures/'
if not os.path.exists(path_fig):
    os.makedirs(path_fig)
filobs = path_out+'observation_MIR' # after input_presimu_*.py

## Read h5 file
spec_unit = read_hdf5(filobs, 'spectral unit')[0]
wvl = read_hdf5(filobs, 'wavelength (microns)')
FnuOBS = read_hdf5(filobs, 'FnuOBS ('+spec_unit+')')
dFnuOBS = read_hdf5(filobs, 'dFnuOBS ('+spec_unit+')')
instr = read_hdf5(filobs, 'spectroscopic module labels')

## Added a NaN at 5 micron
i5 = closest(wvl, 4.3, 'left') + 1
wvl = np.insert(wvl, i5, 4.3)
FnuOBS = np.insert(FnuOBS, i5, np.nan, axis=0)
dFnuOBS = np.insert(dFnuOBS, i5, np.nan, axis=0)

Nw, Ny, Nx = FnuOBS.shape

mode = ['chi2', 'hb']
for m in mode:
    filout = path_out+'presimu_'+m

    ## Read outputs
    ##--------------
    h5_out = pathlib.Path(filout+'.h5')
    if h5_out.exists():
        FnuMOD = read_hdf5(filout, 'FnuMOD ('+spec_unit+')')
        FnuCONT = read_hdf5(filout, 'FnuCONT ('+spec_unit+')')
        FnuLINE = read_hdf5(filout, 'FnuLINE ('+spec_unit+')')
        FnuBAND = read_hdf5(filout, 'FnuBAND ('+spec_unit+')')
        FnuSTAR = read_hdf5(filout, 'FnuSTAR ('+spec_unit+')')
        Pabs = read_hdf5(filout, 'PABS')
        FnuMOD = np.insert(FnuMOD, i5, np.nan, axis=0)
        FnuCONT = np.insert(FnuCONT, i5, np.nan, axis=1)
        FnuLINE = np.insert(FnuLINE, i5, np.nan, axis=1)
        FnuBAND = np.insert(FnuBAND, i5, np.nan, axis=1)
        FnuSTAR = np.insert(FnuSTAR, i5, np.nan, axis=1)
        Pabs = np.insert(Pabs, i5, np.nan, axis=1)

        ## Make plot tables
        Ncont = FnuCONT.shape[0]
        Nline = FnuLINE.shape[0]
        Nband = FnuBAND.shape[0]
        Nstar = FnuSTAR.shape[0]

        ## calib errors
        if m!='chi2':
            ln1pd = read_hdf5(filout, 'Mean of ln(1+delta)') # [Ncalib,Ny,Nx]
            delp1 = np.ones((Nw,Ny,Nx))
            for x in range(Nx):
                for y in range(Ny):
                    delp1[:,y,x] = np.exp(calexpansion(ln1pd[:,y,x], wvl, instr))
            FnuMOD *= delp1
            for i in range(Ncont):
                FnuCONT[i] *= delp1
            for i in range(Nline):
                FnuLINE[i] *= delp1
            for i in range(Nband):
                FnuBAND[i] *= delp1
            for i in range(Nstar):
                FnuSTAR[i] *= delp1
                
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
            Fnu_tab.append((FnuLINE[i,:,:,:]+np.sum(FnuCONT,axis=0)
                            +FnuSTAR[0])*np.prod(Pabs,axis=0))
            if i==0:
                lab_tab.append('line')
            else:
            	lab_tab.append('')
            col_tab.append('g')
        for i in range(Nband):
            Fnu_tab.append((FnuBAND[i,:,:,:]+np.sum(FnuCONT,axis=0)
                            +FnuSTAR[0])*np.prod(Pabs,axis=0))
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
        filename = path_fig+'presimu_'+m
        
        xlab = r'$\lambda\,(\mu m)$'
        if spec_unit=='MKS':
            ylab = r'$F_{\nu}\ \,(W/m2/Hz/sr)$'
        elif spec_unit=='MJyovsr':
            ylab = r'$F_{\nu}\ \,(MJy/sr)$'
        Fnu_tab = np.array(Fnu_tab)
        print(Fnu_tab.shape)
        
        p = pplot(xlog=1, ylog=1,# nonposx='clip', nonposy='clip',
                  xlim=(2.4, 21), ylim=(4e0, 1e3),
                  xlabel=xlab, ylabel=ylab,
                  figsize=(10,8), left=.1, right=.95, top=.95, bottom=.1,
                  titlesize=20, labelsize=20, ticksize=20, title=None, clib=col_tab)
        for i in range(Fnu_tab.shape[0]):
            p.add_plot(wvl, Fnu_tab[i,:,0,0], label=lab_tab[i])
            
        # p = pplot(wvl, FnuOBS/dFnuOBS, 
        #           xlog=0, ylog=0, xlim=(4.5, 22.), #ylim=(0,2.), 
        #           xlab=xlab, ylab=ylab+'/SN ratio', 
        #           legend='best', label='S/N', c='r')
        # p.add_plot(wvl, FnuOBS, label='Data', c='k')
        
        xtic = [2.5, 3, 3.5, 4, 5, 6, 7, 8, 10, 12, 15, 20,]# 30, 40]
        xtic_min = np.arange(2.5, 20., .1)
        p.ax.set_xticks(xtic, minor=False) # major
        p.ax.set_xticks(xtic_min, minor=True) # minor
        p.ax.xaxis.set_major_formatter(ScalarFormatter()) # major
        p.ax.xaxis.set_minor_formatter(NullFormatter()) # minor
        p.ax.legend(loc='upper left', fontsize=20, framealpha=0)
        
        p.save(filename, transparent=True)
        

print('>>> Coucou show_presimu [Done] <<<')

