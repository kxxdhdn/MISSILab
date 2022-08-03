#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""

This is the visualization of fit_[chi2/bb/hb].h5

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
os.makedirs(path_fig, exist_ok=True)
filobs = path_out+'observation_MIR'

## Read h5 file
spec_unit = read_hdf5(filobs, 'spectral unit')[0]
wvl = read_hdf5(filobs, 'wavelength (microns)')
FnuOBS = read_hdf5(filobs, 'FnuOBS ('+spec_unit+')')
dFnuOBS = read_hdf5(filobs, 'dFnuOBS ('+spec_unit+')')
mask = read_hdf5(filobs, 'NaN mask')
instr = read_hdf5(filobs, 'spectroscopic module labels')
Nw, Ny, Nx = FnuOBS.shape
spec_name = read_hdf5(filobs, 'spectrum labels')
# spec_name = np.empty((Ny,Nx), dtype=('<U30'))
# labelist = ['c1_SN10',
#             'c10_SN10',
#             'c1_SN100',
#             'c10_SN100',]
# for y in range(Ny):
#     for x in range(Nx):
#         spec_name[y,x] = labelist[y]+'_'+str(x+1)

## Negative to NaN
# FnuOBS[FnuOBS<0] = np.nan
i5 = closest(wvl, 5, 'left')

mode = ['chi2','bb','hb']
calib = True
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
        ## NaN at 5 micron
        FnuMOD[i5,:,:] = np.nan
        FnuCONT[:,i5,:,:] = np.nan
        FnuLINE[:,i5,:,:] = np.nan
        FnuBAND[:,i5,:,:] = np.nan
        FnuSTAR[:,i5,:,:] = np.nan
        Pabs[:,i5,:,:] = np.nan

        Ncont = FnuCONT.shape[0]
        Nline = FnuLINE.shape[0]
        Nband = FnuBAND.shape[0]
        Nstar = FnuSTAR.shape[0]
        Nextc = Pabs.shape[0]

        ## Mask wvl
        for i in range(Ncont):
            FnuCONT[i][mask] = np.nan
        for i in range(Nline):
            FnuLINE[i][mask] = np.nan
        for i in range(Nband):
            FnuBAND[i][mask] = np.nan
        for i in range(Nstar):
            FnuSTAR[i][mask] = np.nan
        for i in range(Nextc):
            Pabs[i][mask] = np.nan
        Extinction = np.prod(Pabs,axis=0)
        
        ## Spectrum density and calib factors (non-chi2)
        if m!='chi2':
            ## Spectrum density
            n_FnuMOD = read_hdf5(filout, 'n_FnuMOD ('+spec_unit+')')
            ## NaN at 5 micron
            n_FnuMOD[:,i5,:,:] = np.nan
            ## Mask wvl
            for i in range(3):
                n_FnuMOD[i][mask] = np.nan
            
            ## Calib factors
            meanln1pd = read_hdf5(filout, 'Mean of ln(1+delta)') # [Ncalib,Ny,Nx]
            qln1pd = read_hdf5(filout, 'Quantiles of ln(1+delta)')[1] # [Ncalib,Ny,Nx]
            delp1 = np.ones((Nw,Ny,Nx))
            for x in range(Nx):
                for y in range(Ny):
                    delp1[:,y,x] = np.exp(calexpansion(qln1pd[:,y,x], wvl, instr))
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
                
            ## Make plot tables
            ##------------------
            n_Fnu_tab = []
            n_Fnu_tab.append(FnuOBS)
            for i in range(3):
                n_Fnu_tab.append(n_FnuMOD[i,:,:,:])
            n_Fnu_tab = np.array(n_Fnu_tab)

        ## Make plot tables
        ##------------------
        ## Attribute spectra of model (and components) to a table for pplot
        ## Order: obs - total - cont - line - band - star
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
            Fnu_tab.append(FnuCONT[i,:,:,:]*Extinction)
            if i==0:
                lab_tab.append('cont')
            else:
            	lab_tab.append('')
            col_tab.append('pink')
        for i in range(Nline):
            Fnu_tab.append((FnuLINE[i,:,:,:]+np.sum(FnuCONT,axis=0)
                            +FnuSTAR[0])*Extinction)
            if i==0:
                lab_tab.append('line')
            else:
            	lab_tab.append('')
            col_tab.append('g')
        b_Fnu_tab = [] # b
        sumband = 0 # b
        for i in range(Nband):
            Fnu_tab.append((FnuBAND[i,:,:,:]+np.sum(FnuCONT,axis=0)
                            +FnuSTAR[0])*Extinction)
            if i==0:
                lab_tab.append('band')
            else:
                lab_tab.append('')
            col_tab.append('b')
            sumband += FnuBAND[i,:,:,:] # b
        if m!='chi2':
            filmcmc = path_out+'parlog_fit_'+m
            parmcmc = read_hdf5(filmcmc, 'Parameter values')
            t_end = read_hdf5(filmcmc, 'Last index')[0]
            t_burnin = int(t_end/5) - 1
            lnAv = np.median(parmcmc[t_burnin:t_end,-2,:,:],axis=0)
        else:
            lnAv = read_hdf5(filout, 'Best fitted parameter value')[-2,:,:]

        b_Fnu_tab.append( np.exp(np.log(Extinction)/np.exp(lnAv)*1) ) # b
        b_Fnu_tab.append( np.exp(np.log(Extinction)/np.exp(lnAv)*10) ) # b
        b_Fnu_tab.append( np.exp(np.log(Extinction)/np.exp(lnAv)*100) ) # b

        b_Fnu_tab.append(sumband * np.exp(np.log(Extinction)/np.exp(lnAv)*1)) # b
        b_Fnu_tab.append(sumband * np.exp(np.log(Extinction)/np.exp(lnAv)*10)) # b
        b_Fnu_tab.append(sumband * np.exp(np.log(Extinction)/np.exp(lnAv)*100)) # b
        b_Fnu_tab.append(sumband) # b
        for i in range(Nstar):
            Fnu_tab.append(FnuSTAR[i,:,:,:]*Extinction)
            if i==0:
                lab_tab.append('star')
            else:
            	lab_tab.append('')
            col_tab.append('orange')
        
        xlab = r'$\lambda\,(\mu m)$'
        if spec_unit=='MKS':
            ylab = r'$F_{\nu}\ \,(W/m2/Hz/sr)$'
        elif spec_unit=='MJyovsr':
            ylab = r'$F_{\nu}\ \,(MJy/sr)$'
        Fnu_tab = np.array(Fnu_tab)
        b_Fnu_tab = np.array(b_Fnu_tab)
        xtic = [2, 3, 4, 5, 7, 9, 12, 15, 20,]# 30, 40]
        xtic_min = np.arange(2, 20., .1)
        
        ## Plot spectrum density
        ##-----------------------
        if m!='chi2':
            for x in range(Nx):
                for y in range(Ny):
                    if not np.isnan(Fnu_tab[0,:,y,x]).all():
                        
                        p = pplot(#wvl, n_Fnu_tab[0,:,y,x],
                                  xlog=1, ylog=1,# nonposx='clip', nonposy='clip',
                                  ylim=(np.nanmin(abs(Fnu_tab[1,:,y,x]))*.7,
                                        np.nanmax(abs(Fnu_tab[1,:,y,x]))*1.2),
                                  xlabel=xlab, ylabel=ylab, clib=['k','y','r','y'],
                                  xtk=xtic, xtkmi=xtic_min, xtkform='mylog', ytkform='log_sci',
                                  figsize=(10,8), loc='upper left', legendalpha=0,
                                  titlesize=20, xysize=20, tksize=20, legendsize=20)
                        
                        p.ax.fill_between(wvl, n_FnuMOD[0,:,y,x],
                                          n_FnuMOD[2,:,y,x], facecolor='y', alpha=.5)
                        for i in range(3):
                            if i==1:
                                z = 101
                                label = 'Median of models'
                            elif i==0:
                                z = 100
                                label = r'$Median \pm 1-\sigma$'
                            p.add_plot(wvl, n_Fnu_tab[i+1,:,y,x],
                                       zorder=z,
                                       label=label)
                        p.ax.axvline(5.21,c='grey',ls='dashdot',zorder=100) # IRC-SL2
                        p.ax.axvline(7.56,c='grey',ls='dashdot',zorder=100) # SL2-SL1
                        p.ax.axvline(14.28,c='grey',ls='dashdot',zorder=100) # SL1-LL2
                        ## aesthetics
                        for spine in p.ax.spines.items():
                            p.ax.spines[spine[0]].set_linewidth(2) # border width
                        p.ax.xaxis.set_ticks_position('both') # ticks
                        p.ax.yaxis.set_ticks_position('both')
                        p.ax.tick_params(axis='both', which='both', direction='in')
                        p.ax.tick_params(which='major', length=8, width=2)
                        p.ax.tick_params(which='minor', length=4, width=2)

                        filename = path_fig+'fit_'+m+'_density_'+spec_name[y,x]+'.png'
                        p.save(filename, transparent=True, figtight=True)
        
        ## Plot individual fit
        ##---------------------
        for x in range(Nx):
            for y in range(Ny):
                if not np.isnan(Fnu_tab[0,:,y,x]).all():
    
                    p = pplot(xlog=1, ylog=1,# nonposx='clip', nonposy='clip',
                              ylim=(np.nanmin(abs(Fnu_tab[1,:,y,x]))*.7,
                                    np.nanmax(abs(Fnu_tab[1,:,y,x]))*1.2),
                              xlabel=xlab, ylabel=ylab, clib=col_tab,
                              xtk=xtic, xtkmi=xtic_min, xtkform='mylog', ytkform='log_sci',
                              figsize=(16,8), loc='upper left', legendalpha=0,
                              titlesize=20, xysize=20, tksize=20, legendsize=20)
                    for i in range(Fnu_tab.shape[0]):
                        p.add_plot(wvl, Fnu_tab[i,:,y,x], 
                                   zorder=100+i,
                                   label=lab_tab[i])
                    p.ax.axvline(5.21,c='grey',ls='dashdot',zorder=100) # IRC-SL2
                    p.ax.axvline(7.56,c='grey',ls='dashdot',zorder=100) # SL2-SL1
                    p.ax.axvline(14.28,c='grey',ls='dashdot',zorder=100) # SL1-LL2
                    ## Extinction axis
                    ax2 = p.ax.twinx()
                    ax2.plot(wvl, Extinction[:,y,x],
                             c='k', ls='dotted', label=f'Av={}'.format(np.exp(lnAv[:,y,x])))
                    # ax2.set_yscale('log')
                    ax2.set_ylabel('Extinction')
                    # ax2.legend(loc='upper left',fontsize=20,framealpha=0)
                    
                    filename = path_fig+'fit_'+m+'_'+spec_name[y,x]+'.png'
                    p.save(filename, transparent=True, figtight=True)

        
        ## Plot baseline
        ##---------------
        x,y = 0,0
        print(spec_name[y,x])
        
        if not np.isnan(Fnu_tab[0,:,y,x]).all():
            
            p = pplot(xlog=1, ylog=1,# nonposx='clip', nonposy='clip',
                      ylim=(1e-6, 2e0),
                      xlabel=xlab, ylabel=ylab, clib=col_tab,
                      xtk=xtic, xtkmi=xtic_min, xtkform='mylog', ytkform='log_sci',
                      figsize=(16,8), loc='lower center', legendalpha=0,
                      titlesize=20, xysize=20, tksize=20, legendsize=20)
            p.add_plot(wvl, b_Fnu_tab[0,:,y,x], c='k', ls='--', label='Av=1',)
            p.add_plot(wvl, b_Fnu_tab[1,:,y,x], c='k', ls='-.', label='Av=10')
            p.add_plot(wvl, b_Fnu_tab[2,:,y,x], c='k', ls=':', label='Av=100')
            p.add_plot(wvl, b_Fnu_tab[3,:,y,x]/np.nanmax(b_Fnu_tab[6,:,y,x]),
                       c='b', ls='--', label='Bands with Av=1')
            p.add_plot(wvl, b_Fnu_tab[4,:,y,x]/np.nanmax(b_Fnu_tab[6,:,y,x]),
                       c='b', ls='-.', label='Bands with Av=10')
            p.add_plot(wvl, b_Fnu_tab[5,:,y,x]/np.nanmax(b_Fnu_tab[6,:,y,x]),
                       c='b', ls=':', label='Bands with Av=100')
            p.add_plot(wvl, b_Fnu_tab[6,:,y,x]/np.nanmax(b_Fnu_tab[6,:,y,x]),
                       c='b', ls='-', label='Bands without extinction')
            p.ax.tick_params(which='major', length=8, width=2)
            p.ax.tick_params(which='minor', length=4, width=2)

            filename = path_fig+'baseline_'+m+'_'+spec_name[y,x]+'.png'
            p.save(filename, transparent=True, figtight=True)


print('>>> Coucou show_fit [Done] <<<')


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
# p.ax.axvline(cen, c='c')
