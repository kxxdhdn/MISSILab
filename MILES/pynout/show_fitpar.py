#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""

This is the visualization of fitpar (chi2/BB/HB fit)

"""

import os, pathlib
import numpy as np
import matplotlib.pyplot as plt

## astylo
from astylo.mlib import closest
from astylo.iolib import read_hdf5
from astylo.plib import plot2d_m

## local
from utilities import TABand

## Path
mroot = os.path.dirname(os.path.abspath(__file__))+'/../out1/'
path_fig = mroot+'Figures/'
if not os.path.exists(path_fig):
    os.makedirs(path_fig)
filobs = mroot+'galspec'

## Read h5 file
spec_unit = read_hdf5(filobs, 'Spectral unit')[0]
wvl = read_hdf5(filobs, 'Wavelength (microns)')
FnuOBS = read_hdf5(filobs, 'FnuOBS ('+spec_unit+')')
dFnuOBS = read_hdf5(filobs, 'dFnuOBS ('+spec_unit+')')
Ny, Nx = FnuOBS.shape[1:3]

mode = ['chi2','BB','HB']
for m in mode:
    filfit = mroot+'fitpar_'+m

    fitpar = pathlib.Path(filfit+'.h5')
    if fitpar.exists():
        FnuMOD = read_hdf5(filfit, 'FnuMOD ('+spec_unit+')')
        FnuCONT = read_hdf5(filfit, 'FnuCONT ('+spec_unit+')')
        FnuLINE = read_hdf5(filfit, 'FnuLINE ('+spec_unit+')')
        FnuBAND = read_hdf5(filfit, 'FnuBAND ('+spec_unit+')')
        FnuSTAR = read_hdf5(filfit, 'FnuSTAR ('+spec_unit+')')
        Pabs = read_hdf5(filfit, 'PABS')

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
        
        xlab = r'$\lambda\,(\mu m)$'
        if spec_unit=='MKS':
            ylab = r'Surface brightness $\,(W/m2/Hz/sr)$'
        elif spec_unit=='MJyovsr':
            ylab=r'Surface brightness $\,(MJy/sr)$'
        
        Fnu_tab = np.array(Fnu_tab)

        ## Plot individual fit
        ##---------------------
        for x in range(Nx):
            for y in range(Ny):
                title = 'fitpar_'+m+'_('+str(x+1)+','+str(y+1)+')'
                # title = 'fitpar_'+'_('+str(x+1)+','+str(y+1)+')'+m
                
                filename = path_fig+title+'.png'
                                
                p = plot2d_m(wvl, Fnu_tab[:,:,y,x], 
                             xall=wvl, xlog=0, ylog=0, xlim=(4., 22.), #ylim=(0.,1.e3), 
                             title=title, xlab=xlab, ylab=ylab, 
                             lablist=lab_tab, legend='upper left', cl=col_tab)
                
                p.save(filename)

        ## Plot individual band in table
        ##-------------------------------
        for x in range(Nx):
            filename = path_fig+'fitpar_'+m+'_x='+str(x+1)+'.png'
            # filename = path_fig+'x='+str(x+1)+'_fitpar_'+m+'.png'
        
            fig, axes = plt.subplots(nrows=3, ncols=3, figsize=(9,9))
            plt.subplots_adjust(left=.1, bottom=.05, \
                                right=.99, top=.95, wspace=.3, hspace=.4)
            for y in range(Ny):
                px, py = int(y/3), y%3
                axes[px,py].plot(wvl, Fnu_tab[0,:,y,x],c='k')
                for i in range(2):
                    i += 4
                    axes[px,py].plot(wvl, FnuBAND[i,:,y,x],c='b')
                axes[px,py].set_xlim((6., 6.4))
                axes[px,py].set_ylim((0, 6.e2))
                axes[px,py].set_xlabel(xlab)
                axes[px,py].set_ylabel(ylab)
                axes[px,py].set_title(
                    fr'$Spectra ({x+1},{y+1})$')
        
            plt.savefig(filename)
    
# plt.show()

print('>>> Coucou show_fitpar [Done] <<<')

