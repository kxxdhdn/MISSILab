#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""

This is the visualization of fit_[chi2/bb/hb].h5

"""

import os, pathlib
import numpy as np
import matplotlib.pyplot as plt

## rapyuta
from rapyuta.arrays import closest
from rapyuta.inout import read_hdf5
from rapyuta.plots import pplot

## local
from auxil import croot, TABand, calexpansion

## Path
##------
path_out = croot+'/../out/'
path_fig = path_out+'Figures/'
if not os.path.exists(path_fig):
    os.makedirs(path_fig)
filobs = path_out+'observation_MIR'

## Read h5 file
spec_unit = read_hdf5(filobs, 'spectral unit')[0]
spec_name = read_hdf5(filobs, 'spectrum labels')
wvl = read_hdf5(filobs, 'wavelength (microns)')
FnuOBS = read_hdf5(filobs, 'FnuOBS ('+spec_unit+')')
dFnuOBS = read_hdf5(filobs, 'dFnuOBS ('+spec_unit+')')
mask = read_hdf5(filobs, 'NaN mask')
instr = read_hdf5(filobs, 'spectroscopic module labels')
Nw, Ny, Nx = FnuOBS.shape

mode = ['chi2','bb','hb']
for m in mode:
    filout = path_out+'fit_'+m

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
        Ncont = FnuCONT.shape[0]
        Nline = FnuLINE.shape[0]
        Nband = FnuBAND.shape[0]
        Nstar = FnuSTAR.shape[0]
        Nextc = Pabs.shape[0]

        if m!='chi2':
            n_FnuMOD = read_hdf5(filout, 'n_FnuMOD ('+spec_unit+')')
            ## calib errors
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
            for i in range(3):
                n_FnuMOD[i] *= delp1
                
        ## Mask wvl
        for x in range(Nx):
            for y in range(Ny):
                for iw in range(Nw):
                    if mask[iw,y,x]:
                        FnuMOD[iw,y,x] = np.nan
                        for i in range(Ncont):
                            FnuCONT[i,iw,y,x] = np.nan
                        for i in range(Nline):
                            FnuLINE[i,iw,y,x] = np.nan
                        for i in range(Nband):
                            FnuBAND[i,iw,y,x] = np.nan
                        for i in range(Nstar):
                            FnuSTAR[i,iw,y,x] = np.nan
                        for i in range(Nextc):
                            Pabs[i,iw,y,x] = np.nan
                        if m!='chi2':
                            for i in range(3):
                                n_FnuMOD[i,iw,y,x] = np.nan
                
        ## Make plot tables
        ##------------------
        if m!='chi2':
            n_Fnu_tab = []
            n_Fnu_tab.append(FnuOBS)
            for i in range(3):
                n_Fnu_tab.append(n_FnuMOD[i,:,:,:])
            n_Fnu_tab = np.array(n_Fnu_tab)
            
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
        
        xlabel = r'$\lambda\,(\mu m)$'
        if spec_unit=='MKS':
            ylabel = r'Surface brightness $\,(W/m2/Hz/sr)$'
        elif spec_unit=='MJyovsr':
            ylabel=r'Surface brightness $\,(MJy/sr)$'
        
        Fnu_tab = np.array(Fnu_tab)

        ## Plot spectrum density
        ##-----------------------
        if m!='chi2':
            for x in range(Nx):
                for y in range(Ny):
                    if not np.isnan(Fnu_tab[0,:,y,x]).all():
                        # title = 'fit_'+m+'_density_('+str(x+1)+','+str(y+1)+')'
                        title = 'fit_'+m+'_density_'+spec_name[y,x]
                        
                        filename1 = path_fig+title+'.png'
                        
                        p = pplot(wvl, n_Fnu_tab[0,:,y,x],
                                  xlog=1, ylog=1, nonposx='clip', nonposy='clip',
                                  xlim=(2.5, 22.),
                                  ylim=(1.e-2,np.nanmax(Fnu_tab[0,:,y,x])*1.1),
                                  title='M82 '+spec_name[y,x]+' ('+m+', density)',
                                  xlabel=xlabel, ylabel=ylabel,
                                  clib=['k','y','r','y'])
                        for i in range(3):
                            p.add_plot(wvl, n_Fnu_tab[i+1,:,y,x])
                        
                        p.ax.fill_between(wvl, n_FnuMOD[0,:,y,x],
                                          n_FnuMOD[2,:,y,x], facecolor='y')
                        p.save(filename1)
        
        ## Plot individual fit
        ##---------------------
        for x in range(Nx):
            for y in range(Ny):
                if not np.isnan(Fnu_tab[0,:,y,x]).all():
                    # title = 'fit_'+m+'_('+str(x+1)+','+str(y+1)+')'
                    title = 'fit_'+m+'_'+spec_name[y,x]
                    
                    filename2 = path_fig+title+'.png'
    
                    p = pplot(xlog=1, ylog=1, nonposx='clip', nonposy='clip',
                              xlim=(2.5, 22.),
                              ylim=(1.e-2,np.nanmax(Fnu_tab[0,:,y,x])*1.1),
                              title='M82 '+spec_name[y,x]+' ('+m+')', xlabel=xlabel, ylabel=ylabel, 
                              legend='upper left', clib=col_tab)
                    for i in range(Fnu_tab.shape[0]):
                        p.add_plot(wvl, Fnu_tab[i,:,y,x], 
                                   zorder=100+i,
                                   label=lab_tab[i])

                    ## Line cursors
                    # cen = [
                    #       3.3,
                    #       3.45,
                    #       5.2394667,
                    #       5.6437461,
                    #       5.7490305,
                    #       6.0106598, # Small 6.0
                    #       6.2034273,
                    #       6.2672596,
                    #       6.6273788,
                    #       6.8548833,
                    #       7.0791725, # Small 7.1
                    #       7.6000000,
                    #       7.6171123,
                    #       7.8704769,
                    #       8.3623706, # Small 8.3
                    #       8.6204540,
                    #       9.5244838,
                    #       10.707132,
                    #       11.038349,
                    #       11.237893,
                    #       11.400432, # Plateau 11.3
                    #       11.796389,
                    #       11.949674,
                    #       12.626842,
                    #       12.760273,
                    #       13.559342,
                    #       14.257133,
                    #       15.893117,
                    #       16.482868,
                    #       17.082868,
                    #       17.428485,
                    #       17.771096,
                    #       18.925630,
                    #       ]
                    # p.ax.vlines(cen, 0, 1.e4, colors='c')
                    
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
        #         axes[px,py].set_xlabel(xlabel)
        #         axes[px,py].set_ylabel(ylabel)
        #         axes[px,py].set_title(
        #             fr'$Spectra ({x+1},{y+1})$')
        
        #     plt.savefig(filename3)
    
# plt.show()

print('>>> Coucou show_fit [Done] <<<')

