#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""

[presimu] Tracking parameter distribution during Gibbs sampling

"""

import os, pathlib
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter

## rapyuta
from rapyuta.inout import read_hdf5, h5ext
from rapyuta.plots import pplot, plotool

## local
from auxil import croot, TABand

## Path
##------
path_out = croot+'/../out/'
path_fig = path_out+'Figures/'
os.makedirs(path_fig, exist_ok=True)
filobs = path_out+'observation_MIR'
spec_name = read_hdf5(filobs, 'spectrum labels')
instr = read_hdf5(filobs, 'spectroscopic module labels')
chi2_out = pathlib.Path(path_out+'presimu_chi2'+h5ext)
if chi2_out.exists():
    chi2ini = read_hdf5(path_out+'presimu_chi2', 'Best fitted parameter value') # Chi2 result
h5_analysis = path_out+'input_analysis'

## Read h5 file
mode = ['bb']
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
    filmcmc = path_out+'parlog_presimu_'+m

    h5_out = pathlib.Path(filout+h5ext)
    if h5_out.exists():
        # Nmcmc = read_hdf5(filmcmc, 'Length of MCMC')[0]
        Nmcmc = read_hdf5(filmcmc, 'Last index')[0]
        counter = np.arange(Nmcmc)
        # t_end = Nmcmc
        # t_burnin = int(t_end*.3) - 1
        t_end = read_hdf5(h5_analysis, 't_end')[ibm]
        t_burnin = read_hdf5(h5_analysis, 't_burnin')[ibm]
        parname = read_hdf5(filout, 'Parameter label')
        parmcmc = read_hdf5(filmcmc, 'Parameter values')
        
        Npar, Ny, Nx = parmcmc.shape[1:4]
        meanpar = read_hdf5(filout, 'Mean of parameter value')
        # mu = np.mean(parmcmc[t_burnin:t_end,:,:,:], axis=0)
        # print(mu-meanpar)
        stdevpar = read_hdf5(filout, 'Sigma of parameter value')
        # sig = np.std(parmcmc[t_burnin:t_end,:,:,:], axis=0)
        # print(sig-stdevpar)
        # exit()
        qpar = read_hdf5(filout, 'Quantiles of parameter value')

        if calib:
            ln1pdmcmc = read_hdf5(filmcmc, 'ln(1+delta)')
            Ncal = ln1pdmcmc.shape[1]
            meanln1pd = read_hdf5(filout, 'Mean of ln(1+delta)')
            stdevln1pd = read_hdf5(filout, 'Sigma of ln(1+delta)')
            qln1pd = read_hdf5(filout, 'Quantiles of ln(1+delta)')
        
        filmod = path_out+'input_model'
        fixed = read_hdf5(filmod, 'parinfo fixed')
        parini = read_hdf5(filmod, 'parinfo value')
        labB = read_hdf5(filmod, 'label band')
        
        bname = 'Main 11.2'
        ipar = 0
        bex = 0 # Cband: +1 WSband: +2; WLband: +3
        for i, name in enumerate(parname):
            if name[:7]=='lnRband':
                indB = int(name[7:])
                if labB[indB-1]==bname:
                    ipar = i + bex -1

        x,y = 0,0
        print(spec_name[y,x])

            
        ## Parameter distribution (par)
        ##------------------------------
        # for j in range (int(Npar/9)+1):

        #     p = plotool(3, 3, figsize=(12,9))
        #     p.set_fig(wspace=.3, hspace=.4)
        
        #     for jpar in range(9):
        #         jpar += j*9
        #         px, py = int(jpar%9/3), jpar%9%3
        #         if jpar<Npar:
        #             p.set_ax([px,py], xlabel='Parameter value', ylabel='MCMC counts',
        #                      title=fr'${parname[jpar]}\ (BB)\ =\ {qpar[1,jpar,y,x]:.3f}$')
                    
        #             count, bins, patches = p.ax.hist(parmcmc[:Nmcmc,jpar,y,x],
        #                                              bins='sqrt', log=False,
        #                                              facecolor='c', edgecolor='None',
        #                                              alpha=.5,)
        #             ## Chi2 param in red
        #             if fixed[jpar]=='F' and chi2_out.exists():
        #                 p.ax.axvline(chi2ini[jpar,y,x],c='r',label='Chi2',zorder=100)
                        
        #             p.append_handles()
        #             p.set_legend([px,py], loc='upper left', fontsize=10, framealpha=0)
                    
        #     filename = path_fig+'presimu_dist_ipar='+str(j*9+1)+'-'+str(j*9+9)+'.png'
        #     p.save(filename, transparent=True, figtight=True)
        
        
        ## Parameter track (3*3)
        ##-----------------------
        # for j in range(int(Npar/9)+1):
            
        #     p = plotool(3, 3, figsize=(12,9))
        #     p.set_fig(wspace=.3, hspace=.4)
            
        #     for jpar in range(9):
        #         jpar += j*9
        #         px, py = int(jpar%9/3), jpar%9%3
        #         if jpar<Npar:
        #             p.set_ax([px,py], xlabel='Parameter value', ylabel='MCMC counter',
        #                      title=fr'${parname[jpar]}({x+1},{y+1})={qpar[1,jpar,y,x]:.3f}$',)
                    
        #             p.plot(parmcmc[:Nmcmc,jpar,y,x], counter,
        #                    fmt='s', c='c', marker='.', markersize=5,)

        #             ## Fitted (and after burn-in) param in magenta
        #             p.ax.axvline(qpar[1,jpar,y,x],c='m',label=MODE,zorder=100)
        #             p.ax.axvline(qpar[0,jpar,y,x],c='m',ls='dashed',zorder=100)
        #             p.ax.axvline(qpar[2,jpar,y,x],c='m',ls='dashed',zorder=100)
        #             # p.ax.axvline(meanpar[jpar,y,x],c='m',zorder=100)
        #             # p.ax.axvline(meanpar[jpar,y,x]-stdevpar[jpar,y,x],
        #             #            c='m',ls='dashed',zorder=100)
        #             # p.ax.axvline(meanpar[jpar,y,x]+stdevpar[jpar,y,x],
        #             #            c='m',ls='dashed',zorder=100)
        
        #             ## Chi2 param in red
        #             if fixed[jpar]=='F' and chi2_out.exists():
        #                 p.ax.axvline(chi2ini[jpar,y,x],c='r',label='Chi2',zorder=100)

        #             p.append_handles()
        #             p.set_legend([px,py], loc='upper left', fontsize=10, framealpha=0)
            
        #     filename = path_fig+'presimu_ipar='+str(j*9+1)+'-'+str(j*9+9)+'.png'
        #     p.save(filename, transparent=True, figtight=True)


        ## Parameter track (par)
        ##-----------------------
        jj = 0
        for j in range(Npar):

            if fixed[j]=='F':
                jj += 1
                
                p = plotool(2, 1, figsize=(8,8), gridspec_kw={'height_ratios':[3,1]})
                p.set_fig(hspace=.1)

                if chi2_out.exists():
                    xmin = min( np.append(parmcmc[:Nmcmc,j,y,x],chi2ini[j,y,x]) )
                    xmax = max( np.append(parmcmc[:Nmcmc,j,y,x],chi2ini[j,y,x]) )
                else:
                    xmin = min(parmcmc[:Nmcmc,j,y,x])
                    xmax = max(parmcmc[:Nmcmc,j,y,x])
                xgap = xmax - xmin
                xmin += -xgap* .05
                xmax += xgap* .05

                ## walk
                p.set_ax([0,0],xlabel='Parameter value',ylabel='MCMC counter',
                         xlim=(xmin,xmax),
                         xtktoggle=True,xsize=15,ysize=15,xtksize=15,ytksize=15,)
                p.plot(parmcmc[:Nmcmc,j,y,x],counter,
                       fmt='s',c='c',marker='.',markersize=5,label=parname[j])
                ## Chi2 param in red
                if chi2_out.exists():
                    p.ax.axvline(chi2ini[j,y,x],c='r',
                                 label=fr'Chi2: ${chi2ini[j,y,x]:.3f}$',zorder=100)
                ## Fitted (and after burn-in) param in magenta
                p.ax.axvline(qpar[1,j,y,x],c='m',
                             label=MODE+fr': ${qpar[1,j,y,x]:.3f}$',zorder=100)
                p.ax.axvline(qpar[0,j,y,x],c='m',ls='dashed',zorder=100)
                p.ax.axvline(qpar[2,j,y,x],c='m',ls='dashed',zorder=100)
                # p.ax.axvline(meanpar[j,y,x],c='m',zorder=100)
                # p.ax.axvline(meanpar[j,y,x]-stdevpar[j,y,x],c='m',ls='dashed',zorder=100)
                # p.ax.axvline(meanpar[j,y,x]+stdevpar[j,y,x],c='m',ls='dashed',zorder=100)
                p.append_handles()
                p.set_legend(loc='extra upper left',locext=.5,fontsize=15,framealpha=0)
                p.reset_handles()
                # p.ax.xaxis.set_major_formatter(FormatStrFormatter('%.3f'))
                p.ax.invert_yaxis()
                
                ##dist (burn-in)
                p.set_ax([1,0],xlabel='Parameter value',ylabel='MCMC counts',
                         xlim=(xmin,xmax),
                         xsize=15,ysize=15,xtksize=15,ytksize=15,)
                count, bins, patches = p.ax.hist(
                    parmcmc[t_burnin:t_end,j,y,x],bins='sqrt',log=False,
                    facecolor='c',edgecolor='None',alpha=.5,label='After burn-in')
                ## Chi2 param in red
                if chi2_out.exists():
                    p.ax.axvline(chi2ini[j,y,x],c='r',zorder=100)
                ## Fitted (and after burning) param in magenta
                p.ax.axvline(qpar[1,j,y,x],c='m',zorder=100)
                p.ax.axvline(qpar[0,j,y,x],c='m',ls='dashed',zorder=100)
                p.ax.axvline(qpar[2,j,y,x],c='m',ls='dashed',zorder=100)
                p.append_handles()
                p.set_legend(loc='extra upper left',locext=.5,fontsize=15,framealpha=0)

                filename = path_fig+'presimu_par_'+str(jj)+'_'+str(j+1)+'.png'
                p.save(filename, transparent=True, figtight=True)

            
        ## Parameter track (calib)
        ##-------------------------
        if calib:
            for j in range(1,Ncal):
                print(instr[j])
    
                p = plotool(2, 1, figsize=(8,8), gridspec_kw={'height_ratios':[3,1]})
                p.set_fig(hspace=.1)

                xmin = min(np.exp(ln1pdmcmc[:Nmcmc,j,y,x]))
                xmax = max(np.exp(ln1pdmcmc[:Nmcmc,j,y,x]))
                xgap = xmax - xmin
                xmin += -xgap* .05
                xmax += xgap* .05
    
                ## walk
                p.set_ax([0,0],xlabel='Calibration factor',ylabel='MCMC counter',
                         xlim=(xmin,xmax),
                         xtktoggle=True,xsize=15,ysize=15,xtksize=15,ytksize=15,)
                p.plot(np.exp(ln1pdmcmc[:Nmcmc,j,y,x]),counter,
                       fmt='s',c='c',marker='.',markersize=5,label=instr[j],)
                p.ax.axvline(np.exp(qln1pd[1,j,y,x]),c='m',
                             label=MODE+fr': ${np.exp(qln1pd[1,j,y,x]):.3f}$',zorder=100)
                p.ax.axvline(np.exp(qln1pd[0,j,y,x]),c='m',ls='dashed',zorder=100)
                p.ax.axvline(np.exp(qln1pd[2,j,y,x]),c='m',ls='dashed',zorder=100)
                # p.ax.axvline(np.exp(meanln1pd[j,y,x]),c='m',label=MODE,zorder=100)
                # p.ax.axvline(np.exp(meanln1pd[j,y,x]-stdevln1pd[j,y,x]),
                #              c='m',ls='dashed',zorder=100)
                # p.ax.axvline(np.exp(meanln1pd[j,y,x]+stdevln1pd[j,y,x]),
                #              c='m',ls='dashed',zorder=100)
                p.append_handles()
                p.set_legend(loc='extra upper left',locext=.5,fontsize=15,framealpha=0)
                p.reset_handles()
                p.ax.invert_yaxis()
    
                ## dist
                p.set_ax([1,0],xlabel='Calibration factor',ylabel='MCMC counts',
                         xlim=(xmin,xmax),
                         xsize=15, ysize=15, xtksize=15, ytksize=15,)
                count, bins, patches = p.ax.hist(
                    np.exp(ln1pdmcmc[t_burnin:t_end,j,y,x]),bins='sqrt',log=False,
                    facecolor='c',edgecolor='None',alpha=.5,label='After burn-in')
                p.ax.axvline(np.exp(qln1pd[1,j,y,x]),c='m',zorder=100)
                p.ax.axvline(np.exp(qln1pd[0,j,y,x]),c='m',ls='dashed',zorder=100)
                p.ax.axvline(np.exp(qln1pd[2,j,y,x]),c='m',ls='dashed',zorder=100)
                p.append_handles()
                p.set_legend(loc='extra upper left',locext=.5,fontsize=15,framealpha=0)
    
                filename = path_fig+'presimu_calib_'+instr[j]+'.png'
                p.save(filename, transparent=True, figtight=True)

print('>>> Coucou show_presimu_par [Done] <<<')
