#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""

This is the visualization of presimu_[chi2/bb/hb].h5

"""

import os, pathlib
import numpy as np

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
filobs = path_out+'observation_MIR' # after running input_presimu_*.py
h5_analysis = path_out+'input_analysis'

## Read h5 file
spec_unit = read_hdf5(filobs, 'spectral unit')[0]
wvl = read_hdf5(filobs, 'wavelength (microns)')
FnuOBS = read_hdf5(filobs, 'FnuOBS ('+spec_unit+')')
dFnuOBS = read_hdf5(filobs, 'dFnuOBS ('+spec_unit+')')
mask = read_hdf5(filobs, 'NaN mask')
instr = read_hdf5(filobs, 'spectroscopic module labels')

i5 = closest(wvl, 5, 'left')
Nw, Ny, Nx = FnuOBS.shape

## S/N at 15 um
# i15 = closest(wvl, 15)
# print(FnuOBS[i15]/dFnuOBS[i15])
# exit()

mode = ['chi2', 'bb']
calib = True
for m in mode:
    if m=='chi2':
        MODE = 'Chi2'
    elif m=='bb':
        MODE = 'non-HB'
        ibm = 0
    elif m=='hb':
        MODE = 'HB'
        ibm = 1
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

        ## Calib factors (non-chi2)
        if m!='chi2' and calib:
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
            filmcmc = path_out+'parlog_presimu_'+m
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

        ## Plot individual fit
        ##---------------------        
        xlab = r'$\lambda\,(\mu m)$'
        if spec_unit=='MKS':
            ylab = r'$F_{\nu}\ \,(W/m2/Hz/sr)$'
        elif spec_unit=='MJyovsr':
            ylab = r'$F_{\nu}\ \,(MJy/sr)$'
        Fnu_tab = np.array(Fnu_tab)
        xtic = [2, 3, 4, 5, 6, 7, 9, 12, 15, 20,]# 30, 40]
        xtic_min = np.arange(2, 20., .1)
        
        p = pplot(xlog=1, ylog=1,# nonposx='clip', nonposy='clip',
                  ylim=(6e0, 6e3),
                  xlabel=xlab, ylabel=ylab, clib=col_tab,
                  xtk=xtic, xtkmi=xtic_min, xtkform='mylog', ytkform='log_sci',
                  figsize=(10,8), loc='upper left', legendalpha=0,
                  xysize=20, tksize=20, legendsize=20)
        for i in range(Fnu_tab.shape[0]):
            p.add_plot(wvl, Fnu_tab[i,:,0,0], label=lab_tab[i])
            
        # p = pplot(wvl, FnuOBS/dFnuOBS, 
        #           xlog=0, ylog=0, xlim=(4.5, 22.), #ylim=(0,2.), 
        #           xlab=xlab, ylab=ylab+'/SN ratio', 
        #           legend='best', label='S/N', c='r')
        # p.add_plot(wvl, FnuOBS, label='Data', c='k')
        p.ax.axvline(5.35,c='grey',ls='dashdot',zorder=100) # IRC-SL2
        p.ax.axvline(7.56,c='grey',ls='dashdot',zorder=100) # SL2-SL1
        p.ax.axvline(14.28,c='grey',ls='dashdot',zorder=100) # SL1-LL2
        ## Extinction axis
        p.reset_handles()
        if Nextc==1:
            labAv = fr'$A_V={np.exp(lnAv[0,0,0]):.2f}$'
        else:
            labAv = 'Total extinction'
        p.ax = p.ax.twinx()
        p.ax.plot(wvl, Extinction[:,0,0],
                  c='r', ls='dotted', label=labAv)
        p.set_ax(ylog=1,ylim=(1e-2,9e0),ylabel='Extinction',ytkform='mylog',
                 ysize=20,ytksize=20)
        p.ax.legend(loc='upper right',fontsize=20,framealpha=0)
        
        filename = path_fig+'presimu_spec_'+m
        p.save(filename, transparent=True, figtight=True)
        

print('>>> Coucou show_presimu [Done] <<<')

