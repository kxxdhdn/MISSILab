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
h5_analysis = path_out+'input_analysis'

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
plot_den = input("Plot spectrum density (y/n)? ")
for m in mode:
    if m=='chi2':
        MODE = 'Chi2'
    elif m=='bb':
        MODE = 'non-HB'
        ibm = 0
    elif m=='hb':
        MODE = 'HB'
        ibm = 1
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
            if calib:
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
                lab_tab.append('Cont')
            else:
            	lab_tab.append('')
            col_tab.append('pink')
        for i in range(Nline):
            Fnu_tab.append((FnuLINE[i,:,:,:]+np.sum(FnuCONT,axis=0)
                            +FnuSTAR[0])*Extinction)
            if i==0:
                lab_tab.append('Line')
            else:
            	lab_tab.append('')
            col_tab.append('g')
        for i in range(Nband):
            Fnu_tab.append((FnuBAND[i,:,:,:]+np.sum(FnuCONT,axis=0)
                            +FnuSTAR[0])*Extinction)
            if i==0:
                lab_tab.append('Band')
            else:
                lab_tab.append('')
            col_tab.append('b')
        for i in range(Nstar):
            Fnu_tab.append(FnuSTAR[i,:,:,:]*Extinction)
            if i==0:
                lab_tab.append('Star')
            else:
            	lab_tab.append('')
            col_tab.append('orange')

        ## Extinction via lnAv
        parname = read_hdf5(filout, 'Parameter label')
        ind = []
        for i, name in enumerate(parname):
            if name[:4]=='lnAv':
                ind.append(i)
        if m!='chi2':
            filmcmc = path_out+'parlog_fit_'+m
            parmcmc = read_hdf5(filmcmc, 'Parameter values')
            # Nmcmc = read_hdf5(filmcmc, 'Last index')[0]
            # t_end = Nmcmc
            # t_burnin = int(t_end*.3) - 1
            # t_end = read_hdf5(h5_analysis, 't_end')[ibm]
            # t_burnin = read_hdf5(h5_analysis, 't_burnin')[ibm]
            # lnAv = np.median(parmcmc[t_burnin:t_end,ind[0]:ind[-1]+1,:,:],axis=0)
            # lnAv = read_hdf5(filout, 'Mean of parameter value')[ind[0]:ind[-1]+1,:,:]
            lnAv = read_hdf5(filout, 'Quantiles of parameter value')[1,ind[0]:ind[-1]+1,:,:]
        else:
            lnAv = read_hdf5(filout, 'Best fitted parameter value')[ind[0]:ind[-1]+1,:,:]

        ##---------------
        ## Prepare Plots
        ##---------------
        xlab = r'$\lambda\,(\mu m)$'
        if spec_unit=='MKS':
            ylab = r'$F_{\nu}\ \,(W/m2/Hz/sr)$'
        elif spec_unit=='MJyovsr':
            ylab = r'$F_{\nu}\ \,(MJy/sr)$'
        Fnu_tab = np.array(Fnu_tab)
        # xtic = [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 15, 17, 20,]# 30, 40]
        xtic = [2, 3, 4, 5, 6, 7, 9, 12, 15, 20,]# 30, 40]
        xtic_min = np.arange(2, 20, .1)
        
        ## Plot individual fit
        ##---------------------
        for x in range(Nx):
            for y in range(Ny):
                if not np.isnan(Fnu_tab[0,:,y,x]).all():

                    ymin = np.nanmin(abs(Fnu_tab[1,:,y,x])) # obs - 0, mod - 1
                    ymax = np.nanmax(abs(Fnu_tab[1,:,y,x]))
                    # ygap = ymax - ymin
                    # ymin += -ygap* .1
                    # ymax += ygap* .2
                    ymin *= .5
                    ymax *= 5e1
    
                    p = pplot(xlog=1, ylog=1,# nonposx='clip', nonposy='clip',
                              ylim=(ymin,ymax),
                              xlabel=xlab, ylabel=ylab, clib=col_tab,
                              xtk=xtic, xtkmi=xtic_min, xtkform='mylog', ytkform='log_sci',
                              figsize=(10,8), loc='upper left', legendalpha=0,
                              titlesize=20, xysize=20, tksize=20, legendsize=20)
                    for i in range(Fnu_tab.shape[0]):
                        p.add_plot(wvl, Fnu_tab[i,:,y,x], 
                                   zorder=100+i, label=lab_tab[i])
                    p.ax.axvline(5.21,c='grey',ls='dashdot',zorder=100) # IRC-SL2
                    p.ax.axvline(7.56,c='grey',ls='dashdot',zorder=100) # SL2-SL1
                    p.ax.axvline(14.28,c='grey',ls='dashdot',zorder=100) # SL1-LL2
                    ## Extinction axis
                    p.reset_handles()
                    if Nextc==1:
                        labAv = fr'$A_V={np.exp(lnAv[0,y,x]):.2f}$'
                    else:
                        labAv = 'Total extinction'
                    p.ax = p.ax.twinx()
                    p.ax.plot(wvl, Extinction[:,y,x],
                              c='r', ls='dotted', lw='2', label=labAv)
                    p.set_ax(ylog=1, ylabel='Extinction', ytkform='mylog',
                             ylim=(1e-2,5e0),
                             ysize=20, ytksize=20)
                    p.ax.legend(loc='upper right',fontsize=20,framealpha=0)
                    
                    filename = path_fig+'fit_'+m+'_'+spec_name[y,x]+'.png'
                    p.save(filename, transparent=True, figtight=True)
        
        ## Plot spectrum density
        ##-----------------------
        if m!='chi2' and plot_den=='y':
            for x in range(Nx):
                for y in range(Ny):
                    if not np.isnan(Fnu_tab[0,:,y,x]).all():

                        ymin = np.nanmin(abs(n_FnuMOD[0,:,y,x]))
                        ymax = np.nanmax(abs(n_FnuMOD[2,:,y,x]))
                        # ygap = ymax - ymin
                        # ymin += -ygap* .1
                        # ymax += ygap* .2
                        ymin *= .5
                        ymax *= 5e1
                        
                        p = pplot(wvl, Fnu_tab[0,:,y,x], label='Obs',
                                  xlog=1, ylog=1,# nonposx='clip', nonposy='clip',
                                  ylim=(ymin,ymax),
                                  xlabel=xlab, ylabel=ylab, clib=['y','grey','k','grey'],
                                  xtk=xtic, xtkmi=xtic_min, xtkform='mylog', ytkform='log_sci',
                                  figsize=(10,8), loc='upper left', legendalpha=0,
                                  titlesize=20, xysize=20, tksize=20, legendsize=20)
                        p.ax.fill_between(wvl, n_FnuMOD[0,:,y,x], n_FnuMOD[2,:,y,x],
                                          facecolor='grey', alpha=.5, zorder=99)
                        for i in range(3):
                            z = 10
                            if i==1:
                                lw1 = 1
                                z = 100
                                label = 'Median of models'
                            elif i==0:
                                lw1 = 2
                                label = r'Median $\pm 1-\sigma$'
                            else:
                                lw1 = 2
                            p.add_plot(wvl, n_Fnu_tab[i+1,:,y,x],
                                       lw=lw1, zorder=z, label=label)
                        p.ax.axvline(5.21,c='grey',ls='dashdot',zorder=100) # IRC-SL2
                        p.ax.axvline(7.56,c='grey',ls='dashdot',zorder=100) # SL2-SL1
                        p.ax.axvline(14.28,c='grey',ls='dashdot',zorder=100) # SL1-LL2
                        ## Extinction axis
                        p.reset_handles()
                        if Nextc==1:
                            labAv = fr'$A_V={np.exp(lnAv[0,y,x]):.2f}$'
                        else:
                            labAv = 'Total extinction'
                        p.ax = p.ax.twinx()
                        p.ax.plot(wvl, Extinction[:,y,x],
                                  c='r', ls='dotted', lw='2', label=labAv)
                        p.set_ax(ylog=1, ylabel='Extinction', ytkform='mylog',
                                 ylim=(1e-2,5e0),
                                 ysize=20, ytksize=20)
                        p.ax.legend(loc='upper right',fontsize=20,framealpha=0)
                        ## aesthetics
                        # for spine in p.ax.spines.items():
                        #     p.ax.spines[spine[0]].set_linewidth(2) # border width
                        # p.ax.xaxis.set_ticks_position('both') # ticks
                        # p.ax.yaxis.set_ticks_position('both')
                        # p.ax.tick_params(axis='both', which='both', direction='in')
                        # p.ax.tick_params(which='major', length=8, width=2)
                        # p.ax.tick_params(which='minor', length=4, width=2)

                        filename = path_fig+'fit_'+m+'_density_'+spec_name[y,x]+'.png'
                        p.save(filename, transparent=True, figtight=True)


print('>>> Coucou show_fit [Done] <<<')
