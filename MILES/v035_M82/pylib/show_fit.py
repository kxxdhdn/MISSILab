#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""

This is the visualization of fit_[chi2/bb/hb].h5

"""

import os, pathlib
import numpy as np
import matplotlib.pyplot as plt

## laputan
from laputan.arrays import closest
from laputan.inout import read_hdf5
from laputan.plots import pplot

## local
from librarian import croot, TABand

## Path
##------
path_out = croot+'/../out/'
path_fig = path_out+'Figures/'
if not os.path.exists(path_fig):
    os.makedirs(path_fig)
filobs = path_out+'observation_MIR' # after input_[src]_[chi2/bb/hb].py

## Read h5 file
spec_unit = read_hdf5(filobs, 'spectral unit')[0]
wvl = read_hdf5(filobs, 'wavelength (microns)')
FnuOBS = read_hdf5(filobs, 'FnuOBS ('+spec_unit+')')
dFnuOBS = read_hdf5(filobs, 'dFnuOBS ('+spec_unit+')')
Ny, Nx = FnuOBS.shape[1:3]

labelist = ['c1_SN5','c3_SN5','c9_SN5',
            'c1_SN20','c3_SN20','c9_SN20',
            'c1_SN80','c3_SN80','c9_SN80',]

mode = ['chi2','bb','hb']
for m in mode:
    filout = path_out+'fit_'+m

    h5_out = pathlib.Path(filout+'.h5')
    if h5_out.exists():
        FnuMOD = read_hdf5(filout, 'FnuMOD ('+spec_unit+')')
        FnuCONT = read_hdf5(filout, 'FnuCONT ('+spec_unit+')')
        FnuLINE = read_hdf5(filout, 'FnuLINE ('+spec_unit+')')
        FnuBAND = read_hdf5(filout, 'FnuBAND ('+spec_unit+')')
        FnuSTAR = read_hdf5(filout, 'FnuSTAR ('+spec_unit+')')
        Pabs = read_hdf5(filout, 'PABS')
        if m!='chi2':
            n_FnuMOD = read_hdf5(filout, 'n_FnuMOD ('+spec_unit+')')

        ## Make plot tables
        if m!='chi2':
            n_Fnu_tab = []
            n_Fnu_tab.append(FnuOBS)
            for i in range(3):
                n_Fnu_tab.append(n_FnuMOD[i,:,:,:])
            n_Fnu_tab = np.array(n_Fnu_tab)
            
        Ncont = FnuCONT.shape[0]
        Nline = FnuLINE.shape[0]
        Nband = FnuBAND.shape[0]
        Nstar = FnuSTAR.shape[0]
        Fnu_tab = []
        lab_tab = []
        col_tab = ['k']
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

        ## Plot spectrum density
        ##-----------------------
        if m!='chi2':
            for x in range(Nx):
                for y in range(Ny):
                    title = 'fit_'+m+'_density_('+str(x+1)+','+str(y+1)+')'
                    # title = 'fit_'+'_('+str(x+1)+','+str(y+1)+')'+m
                    
                    filename1 = path_fig+title+'.png'

                    p = pplot(wvl, n_Fnu_tab[0,:,y,x],
                              xlog=0, ylog=0, xlim=(4., 22.), #ylim=(0.,1.e3),
                              title=title+' -['+labelist[y]+']', xlab=xlab, ylab=ylab,
                              clib=['k','y','r','y'])
                    for i in range(3):
                        p.add_plot(wvl, n_Fnu_tab[i,:,y,x])
    
                    p.ax.fill_between(wvl, n_FnuMOD[0,:,y,x], n_FnuMOD[2,:,y,x], facecolor='y')
                    p.save(filename1)
        
        ## Plot individual fit
        ##---------------------
        for x in range(Nx):
            for y in range(Ny):
                title = 'fit_'+m+'_('+str(x+1)+','+str(y+1)+')'
                # title = 'fit_'+'_('+str(x+1)+','+str(y+1)+')'+m
                
                filename2 = path_fig+title+'.png'

                p = pplot(xlog=0, ylog=0, xlim=(4., 22.), #ylim=(0.,1.e3), 
                          title=title+' -['+labelist[y]+']', xlab=xlab, ylab=ylab, 
                          legend='upper left', clib=col_tab)
                for i in range(Fnu_tab.shape[0]):
                    p.add_plot(wvl, Fnu_tab[i,:,y,x], label=lab_tab[i])
                
                p.save(filename2)
                
        ## Plot individual band in table
        ##-------------------------------
        # for x in range(Nx):
        #     filename3 = path_fig+'fit_'+m+'_x='+str(x+1)+'.png'
        #     # filename3 = path_fig+'x='+str(x+1)+'_fit_'+m+'.png'
        
        #     fig, axes = plt.subplots(nrows=3, ncols=3, figsize=(9,9))
        #     plt.subplots_adjust(left=.1, bottom=.05, \
        #                         right=.99, top=.95, wspace=.3, hspace=.4)
        #     for y in range(Ny):
        #         px, py = int(y/3), y%3
        #         axes[px,py].plot(wvl, Fnu_tab[0,:,y,x],c='k')
        #         for i in range(2):
        #             i += 4
        #             axes[px,py].plot(wvl, FnuBAND[i,:,y,x],c='b')
        #         axes[px,py].set_xlim((6., 6.4))
        #         axes[px,py].set_ylim((0, 6.e2))
        #         axes[px,py].set_xlabel(xlab)
        #         axes[px,py].set_ylabel(ylab)
        #         axes[px,py].set_title(
        #             fr'$Spectra ({x+1},{y+1})$')
        
        #     plt.savefig(filename3)
    
# plt.show()

print('>>> Coucou show_fit [Done] <<<')

